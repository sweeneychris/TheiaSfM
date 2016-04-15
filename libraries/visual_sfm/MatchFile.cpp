////////////////////////////////////////////////////////////////////////////
//	File:		MatchFile.cpp
//	Author:		Changchang Wu (ccwu1130@gmail.com)
//	Description :
//
//  Copyright (c) 2011  Changchang Wu (ccwu1130@gmail.com)
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation; either
//  Version 3 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////

#include "stdio.h"
#include <cassert>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <iomanip>
using namespace std;
#include "MatchFile.h"

#ifndef O_BINARY
#define O_BINARY  0
#define O_TEXT    0
#endif

#ifndef DISABLE_NAMED_MUTEX
#include "NamedMutex.h"
#endif

namespace visual_sfm {

//////////////////////////////////////
int MatchFile::multi_thread_mode = 0;
int MatchFile::file_title_mode = 0;
int MatchFile::num_match_verification = 1;
int MatchFile::delay_header_update = 0;
int MatchFile::record_reservation = 20;

// function for printing out messages.
int (*MatchFile::matchfile_printf)(const char* format, ...) = ::printf;
#define printf matchfile_printf

#ifdef _WIN64
#define lseek _lseeki64
#define strcasecmp _stricmp
#endif

#ifdef WIN32
#define PATH_SLASH '\\'
#define PATH_SLASH_X '/'
#define PATH_PARENT "..\\"
#else
#define PATH_SLASH '/'
#define PATH_SLASH_X '\\'
#define PATH_PARENT "../"
#endif

void MatchFile::SetPrintFunction(int (*printf_func)(const char* format, ...)) {
  matchfile_printf = printf_func;
}

MatchFile::MatchFile() {
  _fid = 0;
  _header.definition_size = 0;
  _header.file_count = 0;
  _header.feature_count = 0;
  _recordloc_unchanged = 0;
  _header_changed = 0;
  _opened_forwrite = 0;
  _mutex = NULL;
  _matchfile_mode = 0;
  _title = _filepath;
}

MatchFile::~MatchFile() { CloseMatchFile(); }

MatchFile::MatchFile(const char* image_path, int mode) {
  _fid = 0;
  _mutex = NULL;
  _header.definition_size = 0;
  _header.file_count = 0;
  _header.feature_count = 0;
  _recordloc_unchanged = 0;
  _header_changed = 0;
  OpenMatchFile(image_path, 0, mode);
}

int MatchFile::IsMatchFileOfImage(const char* fullpath, const char* title) {
  if (file_title_mode && title)
    return strcmp(_title, title) == 0;
  else
    return strcmp(_filepath, fullpath) == 0;
}

void MatchFile::GetMatchFilePath(char* match_path) {
  strcpy(match_path, _filepath);

  if (_matchfile_mode == 0)
    strcat(match_path, ".mat");
  else if (_matchfile_mode == 1)
    strcat(match_path, ".gmat");
  else
    strcat(match_path, ".xmat");

  //////////////////////
  _title = strrchr(_filepath, PATH_SLASH);
  if (_title == NULL)
    _title = _filepath;
  else
    _title++;
}

void MatchFile::GetMatchFolderPath(char* match_path) {
  strcpy(match_path, _filepath);
  char* p = strrchr(match_path, PATH_SLASH);
  p[1] = 0;
}

int MatchFile::MakeWritable() {
  char match_path[MAX_PATH];
  // already opened for write
  if (_fid <= 0) return 0;
  if (_opened_forwrite) return 1;
  //////////////////////////////

  GetMatchFilePath(match_path);

  close(_fid);
  _fid =
      open(match_path, O_BINARY | O_CREAT | O_RDWR, S_IWRITE | S_IREAD);
  if (_fid > 0) {
    _opened_forwrite = 1;
    return 1;
  } else {
    _opened_forwrite = 0;
    return 0;
  }
}

void MatchFile::VerifyFeatureCount(int feature_count) {
  if (_fid <= 0) return;

  if (_header.feature_count <= 0) {
    // update the feature count
    if (feature_count > 0) {
      _header.feature_count = feature_count;
      if (MakeWritable()) {
        lseek(_fid, 0, SEEK_SET);
        write(_fid, &_header, sizeof(_header));
        // printf("Update feature count [%d]!\r\n", feature_count);
      }
    }
  } else if (feature_count > 0 && _header.feature_count != feature_count) {
    // destroy all match record.
    DeleteMatchFile();
    printf("Delete match record: feature changed!\r\n");
  }
}

void MatchFile::DeleteMatchRecord(const char* pattern) {
  if (_fid <= 0) return;
  int changed = 0;
  for (int i = 0; i < _match_records.size(); ++i) {
    if (strstr(_match_records[i]->file_name, pattern)) {
      int NM;
      RecordLoc* loc = _match_records[i];
      lseek(_fid, loc->read_loc, SEEK_SET);
      read(_fid, &NM, sizeof(int));

      ////////////////////////////////
      if (NM == 0) continue;
      MakeWritable();
      lseek(_fid, loc->read_loc, SEEK_SET);
      NM = 0;
      write(_fid, &NM, sizeof(int));
      changed = 1;
    }
  }
}

void MatchFile::DeleteMatchFile() {
  char match_path[MAX_PATH];
  GetMatchFilePath(match_path);
  CloseMatchFile();
  remove(match_path);
}

int MatchFile::OpenMatchFile(const char* imagepath, int write, int mode) {
  char match_path[MAX_PATH];

  CloseMatchFile();

  //
  _matchfile_mode = mode;

  ///////////////
  strcpy(_filepath, imagepath);

  GetMatchFilePath(match_path);

#ifndef DISABLE_NAMED_MUTEX
  if (multi_thread_mode) {
    // use file title as the mutex name
    const char* slash = strrchr(match_path, PATH_SLASH);
    _mutex = new NamedMutex(slash ? slash + 1 : match_path);
  } else {
    _mutex = NULL;
  }
#endif

  if (!write)
    _fid = open(match_path, O_BINARY | O_RDONLY);
  else
    _fid =
        open(match_path, O_BINARY | O_CREAT | O_RDWR, S_IWRITE | S_IREAD);

  if (_fid > 0) {
    // read file _header
    int readsz = read(_fid, (char*)(&_header), sizeof(_header));
    _opened_forwrite = write;

    if (sizeof(_header) == readsz && _header.definition_size > 0 &&
        _header.definition_buf >= _header.definition_size &&
        _header.version == MATCH_FILE_VERSION_3) {
      assert(_header.file_count >= 0);
      _record_definition.resize(_header.definition_buf, 0);
      read(_fid, (char*)&_record_definition[0], _header.definition_size);
      GetRecordList();
      _header_changed = 0;
      _recordloc_unchanged = _header.definition_size;
    } else if (readsz > 0) {
      printf("Unsupported MatchFile [%s]\r\n", _filepath);
      ResetMatchFile();
      if (write) UpdateHeaderAndRecordLoc();
    }
    return 1;
  } else {
#ifndef DISABLE_NAMED_MUTEX
    if (_mutex) {
      delete _mutex;
      _mutex = NULL;
    }
#endif
    _filepath[0] = 0;
    return 0;
  }
}

void MatchFile::GetRecordList() {
  char* p = &_record_definition[0];
  _match_records.resize(_header.file_count);
  for (int i = 0; i < _match_records.size(); ++i) {
    RecordLoc* rp = (RecordLoc*)p;
    _match_records[i] = rp;

    if (rp->read_loc <= 0 || rp->block_size == 0 || rp->extra_size > MAX_PATH ||
        (rp->fcount > 1000000 || rp->fcount < 0)) {
      _match_records.resize(i);
      _header.file_count = i;
      _header.definition_size = (p - &_record_definition[0]);
      break;
    } else {
      p += sizeof(RecordLoc) + _match_records[i]->extra_size;
      // convert file name for different operation system
      for (char* c = rp->file_name; *c; ++c) {
        if (*c == PATH_SLASH_X) *c = PATH_SLASH;
      }
    }
  }
}

void MatchFile::GetMatchedImageList(vector<string>& paths) {
  char match_path[MAX_PATH];
  GetMatchFolderPath(match_path);
  string folder = match_path;
  paths.resize(0);

  for (int i = 0; i < _match_records.size(); ++i) {
    paths.push_back(folder + _match_records[i]->file_name);
  }
}

void MatchFile::SaveSubsetMatch(vector<string>& paths, vector<int>& fc) {
  MatchFile sub;
  sub.OpenMatchFile(_filepath, 1, 2);
  sub.ResetMatchFile();
  ////
  int NM;
  Points<int> matches;
  TwoViewGeometry tvg;
  Points<int> inliers;

  /////////////////////////////
  for (int i = 0; i < paths.size(); ++i) {

    if (GetPMatch(paths[i].c_str(), fc[i], NM, matches)) {
      sub.WritePMatch(paths[i].c_str(), fc[i], NM, matches, 0);
      if (GetIMatch(paths[i].c_str(), fc[i], tvg, inliers)) {
        sub.WriteIMatch(paths[i].c_str(), 0, NM, tvg, inliers);
      }
    }
  }
}

int MatchFile::HaveMatchRecord(const char* absolute_path) {
  char relative_path[MAX_PATH];
  GetRelativePath(absolute_path, relative_path);
  return GetImageIndex(relative_path) >= 0;
}

int MatchFile::GetImageIndex(const char* relative_path) {
  if (file_title_mode) {
    const char* ps = strrchr(relative_path, PATH_SLASH);
    ps = (ps == NULL ? relative_path : ps + 1);

    for (int i = 0; i < _match_records.size(); i++) {
      RecordLoc* loc = _match_records[i];
      const char* p = strrchr(loc->file_name, PATH_SLASH);
      p = (p == NULL) ? loc->file_name : p + 1;
      if (strcmp(p, ps) == 0) return i;
    }
  } else {
    for (int i = 0; i < _match_records.size(); i++) {
      if (strcmp(_match_records[i]->file_name, relative_path) == 0) return i;
    }
  }
  return -1;
}

void MatchFile::GetRelativePath(const char* path1, const char* path2,
                                char* relative_path) {

  int i = 0;
  while (path2[i] == path1[i] && path2[i] && path1[i]) i++;
  while (i > 0 && path1[i] != PATH_SLASH) i--;
  if (i == 0) {
    strcpy(relative_path, path2);
  } else {
    relative_path[0] = 0;
    const char* p = path1 + i + 1;
    char* rp = relative_path;
    int rc1 = 0, rc2 = 0;
    while (*p) {
      if (p[0] == PATH_SLASH && p[1] != PATH_SLASH) {
        strcpy(rp, PATH_PARENT);
        rp += 3;
        rc1++;
      }
      p++;
    }
    p = path2 + i + 1;
    while (*p) {
      if (p[0] == PATH_SLASH && p[1] != PATH_SLASH) rc2++;
      p++;
    }
    if (rc1 <= 1 || rc2 <= 1)
      strcpy(rp, path2 + i + 1);
    else
      strcpy(relative_path, path2);
  }
}

inline void MatchFile::GetRelativePath(const char* image_match,
                                       char* relative_path) {
  if (file_title_mode) {
    const char* pslash = strrchr(image_match, PATH_SLASH);
    if (pslash)
      strcpy(relative_path, pslash + 1);
    else
      strcpy(relative_path, image_match);
  } else {
    GetRelativePath(_filepath, image_match, relative_path);
  }
}

void MatchFile::MoveRecordToEnd(int i) {
  if (i >= _header.file_count - 1 || i < 0 || _header.file_count < 2) return;
  char* p0 = (char*)_match_records[0];
  char* p1 = (char*)_match_records[i];
  char* p2 = (char*)_match_records[i + 1];
  size_t loc1 = p1 - p0;
  size_t loc2 = p2 - p0;
  int len = loc2 - loc1;
  vector<char> temp;
  temp.insert(temp.begin(), p1, p2);
  _record_definition.erase(_record_definition.begin() + loc1,
                           _record_definition.begin() + loc2);
  //_record_definition.insert(_record_definition.end(), temp.begin(),
  //temp.begin() + len);
  _record_definition.insert(
      _record_definition.begin() + _header.definition_size - len, temp.begin(),
      temp.begin() + len);
  _recordloc_unchanged = min(_recordloc_unchanged, loc1);
  GetRecordList();

  if (i == 0) {
    _header.definition_buf += len;
    _header_changed = 1;
  } else {
    _match_records[i - 1]->trash_size += len;
  }
  //

  _match_records[_header.file_count - 1]->read_loc =
      _match_records[_header.file_count - 2]->read_loc +
      _match_records[_header.file_count - 2]->block_size;

  printf("MoveRecordToEnd [%s]\r\n", _filepath);
}

int MatchFile::AddImageMatch(const char* relative_path) {
  // see if the name buffer has enough space..
  int slen = (strlen(relative_path) / 4 + 1) * 4;
  int elen = slen - 4;
  int rlen = elen + sizeof(RecordLoc);

  //////////////////
  _header_changed = 1;
  //////////////
  if (_header.file_count == 0) {
    _header.version = MATCH_FILE_LATEST_VERSION;
    _header.file_count = 1;
    _header.definition_buf =
        max(20, record_reservation) * (MAX_PATH + sizeof(RecordLoc));
    _record_definition.resize(_header.definition_buf, 0);
    _header.definition_size = rlen;
    RecordLoc* loc = (RecordLoc*)(&_record_definition[0]);
    _match_records.resize(1);
    _match_records[0] = loc;
    strcpy(loc->file_name, relative_path);
    loc->extra_size = elen;
    loc->trash_size = 0;
    loc->block_size = 0;
    loc->read_loc = sizeof(FileHeader) + _header.definition_buf;
    // should write the _header now?
    return 0;
  } else if (_header.definition_size + rlen < _header.definition_buf) {
    RecordLoc* rec = (RecordLoc*)(&_record_definition[_header.definition_size]);
    rec->block_size = 0;
    rec->extra_size = elen;
    rec->read_loc =
        _match_records.back()->read_loc + _match_records.back()->block_size;
    strcpy(rec->file_name, relative_path);
    _match_records.push_back(rec);
    _header.definition_size += rlen;
    _header.file_count++;
    return _header.file_count - 1;
  } else {
    // size_t szr = max(((int) (_match_records.size() + 1)), record_reservation)
    // * (MAX_PATH + sizeof(RecordLoc));
    int extra_needed = rlen + _header.definition_buf;
    int match_move_needed = 0;
    int extra_space = 0;

    // move the match data
    while (extra_space < extra_needed) {
      RecordLoc* loc = _match_records[match_move_needed];
      extra_space += (loc->block_size + loc->trash_size);
      match_move_needed++;
    }
    assert(match_move_needed > 0);

    // read match data that needs to move
    vector<char> read_buffer(extra_space);
    char* pbuf = &read_buffer[0];
    lseek(_fid, _match_records[0]->read_loc, SEEK_SET);
    read(_fid, &read_buffer[0], extra_space);

    // write match data to the end of file
    int matchwrite_location =
        _match_records.back()->read_loc + _match_records.back()->block_size;
    lseek(_fid, matchwrite_location, SEEK_SET);
    for (int i = 0; i < match_move_needed; i++) {
      RecordLoc* loc = _match_records[i];
      write(_fid, pbuf, loc->block_size);
      pbuf += (loc->block_size + loc->trash_size);
      loc->read_loc = matchwrite_location;
      matchwrite_location += loc->block_size;
      loc->trash_size = 0;
    }

    // organize file names
    int definition_move =
        ((char*)_match_records[match_move_needed]) - &_record_definition[0];
    if (match_move_needed < _header.file_count) {
      // reorganize strings..move several string to the end
      vector<char> temp = _record_definition;
      std::copy(temp.begin() + definition_move,
                temp.begin() + _header.definition_size,
                _record_definition.begin());
      std::copy(temp.begin(), temp.begin() + definition_move,
                _record_definition.begin() + _header.definition_size -
                    definition_move);
    }
    _header.definition_buf = _header.definition_buf * 2 + rlen;
    _record_definition.resize(_header.definition_buf, 0);
    GetRecordList();

    RecordLoc* rec = (RecordLoc*)(&_record_definition[_header.definition_size]);
    rec->read_loc = matchwrite_location;
    rec->trash_size = 0;
    rec->block_size = 0;
    rec->extra_size = elen;
    _match_records.push_back(rec);
    strcpy(rec->file_name, relative_path);

    // write all previous data to file
    lseek(_fid, sizeof(_header), SEEK_SET);
    write(_fid, &_record_definition[0], _header.definition_size);

    ////
    _header.definition_size += rlen;
    _header.file_count++;
    _recordloc_unchanged = 0;
    return _header.file_count - 1;
  }
}

void MatchFile::ResetMatchFile() {
  _header.file_count = 0;
  _header.definition_size = 0;
  _header.definition_buf = 0;
  _header.version = MATCH_FILE_LATEST_VERSION;
  _match_records.resize(0);
  _record_definition.resize(0);
  _recordloc_unchanged = 0;
  _header_changed = 1;
}

void MatchFile::CloseMatchFile() {
  if (_fid > 0) {
    if (_opened_forwrite && delay_header_update) UpdateHeaderAndRecordLoc();
    ResetMatchFile();
    close(_fid);
    _fid = 0;
    _opened_forwrite = 0;
    _matchfile_mode = 0;
  }
#ifndef DISABLE_NAMED_MUTEX
  if (_mutex) {
    delete _mutex;
    _mutex = NULL;
  }
#endif
}

void MatchFile::WriteFMatch(const char* image1, const char* image2, int NF,
                            int* index1, int* index2, float F[3][3]) {
  if (index1 == NULL || index2 == NULL) return;
  int* indices[2] = {index1, index2};
  Points<int> inliers(indices, NF, 2);
  TwoViewGeometry tvg(NF, F);
  WriteIMatch(image1, image2, 0, tvg, inliers);
}

void MatchFile::WriteIMatch(MatchFile* mat, const char* image1,
                            const char* image2, int NM, TwoViewGeometry& tvg,
                            Points<int>& inliers) {

  if (mat == NULL || !mat->IsValid()) {
    WriteIMatch(image1, image2, NM, tvg, inliers);
  } else if (strcasecmp(mat->_filepath, image1) == 0) {
    if (mat->MakeWritable()) {
      mat->WriteIMatch(image2, 0, NM, tvg, inliers);
    }
    MatchFile file2;
    if (file2.OpenMatchFile(image2, 1)) {
      file2.WriteIMatch(image1, 1, NM, tvg, inliers);
    }
  } else if (strcasecmp(mat->_filepath, image2) == 0) {
    if (mat->MakeWritable()) {
      mat->WriteIMatch(image1, 1, NM, tvg, inliers);
    }

    MatchFile file1;
    if (file1.OpenMatchFile(image1, 1)) {
      file1.WriteIMatch(image2, 0, NM, tvg, inliers);
    }
  }
}

void MatchFile::WriteIMatch(const char* image1, const char* image2, int NM,
                            TwoViewGeometry& tvg, Points<int>& inliers) {

  MatchFile file1, file2;
  if (file1.OpenMatchFile(image1, 1)) {
    file1.WriteIMatch(image2, 0, NM, tvg, inliers);
    file1.CloseMatchFile();
  }

  if (file2.OpenMatchFile(image2, 1)) {
    file2.WriteIMatch(image1, 1, NM, tvg, inliers);
    file2.CloseMatchFile();
  }
}

void MatchFile::WriteIMatch(const char* image_match, int reverse, int NM,
                            TwoViewGeometry& tvg, Points<int>& inliers) {
  int index, fnm, szmin;
  char relative_path[MAX_PATH];
  if (image_match == NULL || image_match[0] == 0) return;
  GetRelativePath(image_match, relative_path);
  index = GetImageIndex(relative_path);
  if (index < 0) return;

  RecordLoc* loc = _match_records[index];
  MatchRecordV3 rec;

  assert(loc->read_loc > 0);

  // read in the number of
  if (num_match_verification || NM == 0) {
    lseek(_fid, loc->read_loc, SEEK_SET);
    read(_fid, &fnm, sizeof(int));
    NM = fnm;
  } else {
    lseek(_fid, loc->read_loc + sizeof(int), SEEK_SET);
  }

  assert(tvg.NF <= 0 || tvg.NE == 0 || tvg.NF == tvg.NE);

  szmin = (1 + (NM + tvg.NF) * 2) * sizeof(int) + sizeof(rec);

  if (szmin > loc->block_size)  // this is unlikey to happen, right?
  {
    // not enough size.
    size_t offset = ((char*)loc) - ((char*)_match_records[0]);
    _recordloc_unchanged = min(_recordloc_unchanged, offset);
    if (loc->block_size == 0 || index == _header.file_count - 1) {
      int bs = sizeof(MatchRecordV3) + max(NM + tvg.NF, 8) * sizeof(int) * 4;
      // last one in the file
      loc->block_size = bs;
    } else if (szmin > loc->block_size + loc->trash_size) {
      // not enough space?
      MoveRecordToEnd(index);
      index = _header.file_count - 1;
      loc = _match_records[index];
      loc->block_size =
          sizeof(MatchRecordV3) + max(NM + tvg.NF, 8) * sizeof(int) * 4;
    } else {
      loc->block_size += loc->trash_size;
      loc->trash_size = 0;
    }

    if (!delay_header_update) {
      printf("Reallocate space fo record %d\r\n", index);
      UpdateHeaderAndRecordLoc();
      lseek(_fid, loc->read_loc + sizeof(int), SEEK_SET);
    }
  }

  //
  // Write inlier pairs
  if (reverse == 0) {
    rec.tvg.SetGeometry(tvg);
    // write inlier match record
    write(_fid, &rec, sizeof(rec));
    lseek(_fid, NM * sizeof(int) * 2, SEEK_CUR);
    if (tvg.NF > 0) {
      write(_fid, inliers[0], sizeof(int) * tvg.NF);
      write(_fid, inliers[1], sizeof(int) * tvg.NF);
    }

  } else {
    rec.tvg.SetGeometryR(tvg);
    ;
    // write inlier match record
    write(_fid, &rec, sizeof(rec));

    lseek(_fid, NM * sizeof(int) * 2, SEEK_CUR);
    if (tvg.NF > 0) {
      write(_fid, inliers[1], sizeof(int) * tvg.NF);
      write(_fid, inliers[0], sizeof(int) * tvg.NF);
    }
  }
}

void MatchFile::WritePMatch(const char* image1, const char* image2, int fc1,
                            int fc2, int NM, Points<int>& matches) {
  MatchFile file1, file2;
  if (file1.OpenMatchFile(image1, 1)) {
    file1.WritePMatch(image2, fc2, NM, matches, 0);
    file1.CloseMatchFile();
  }

  if (file2.OpenMatchFile(image2, 1)) {
    file2.WritePMatch(image1, fc1, NM, matches, 1);
    file2.CloseMatchFile();
  }
}

void MatchFile::WritePMatch(const char* image_match, int FC, int NM,
                            Points<int>& matches, int reverse) {
  int index;
  char relative_path[MAX_PATH];
  if (image_match == NULL || image_match[0] == 0) return;
  GetRelativePath(image_match, relative_path);
  index = GetImageIndex(relative_path);
  if (index < 0) index = AddImageMatch(relative_path);
  RecordLoc* loc = _match_records[index];
  loc->fcount = FC;
  assert(loc->read_loc > 0);
  size_t offset = ((char*)loc) - ((char*)_match_records[0]);
  _recordloc_unchanged = min(_recordloc_unchanged, offset);
  if (loc->block_size == 0 || index == _header.file_count - 1) {
    int bs = sizeof(MatchRecordV3) + max(NM, 8) * sizeof(int) * 6;
    // last one in the file
    loc->block_size = bs;

  } else {
    int bsmin = sizeof(MatchRecordV3) + max(NM, 0) * sizeof(int) * 4;
    if (bsmin > loc->block_size) {
      if (bsmin > loc->block_size + loc->trash_size) {
        // not enough space?
        MoveRecordToEnd(index);
        index = _header.file_count - 1;
        loc = _match_records[index];
        loc->block_size = sizeof(MatchRecordV3) + max(NM, 8) * sizeof(int) * 6;
      } else {
        loc->block_size += loc->trash_size;
        loc->trash_size = 0;
      }
    }
  }

  // update file _header and locs
  if (!delay_header_update) {
    UpdateHeaderAndRecordLoc();
  }

  MatchRecordV3 rec;

  // write pmatch
  lseek(_fid, loc->read_loc, SEEK_SET);
  // write match record
  write(_fid, &NM, sizeof(int));
  write(_fid, &rec, sizeof(rec));

  // write matched pairs
  if (NM <= 0) return;
  if (reverse == 0) {
    write(_fid, matches[0], sizeof(int) * NM);
    write(_fid, matches[1], sizeof(int) * NM);
  } else {
    write(_fid, matches[1], sizeof(int) * NM);
    write(_fid, matches[0], sizeof(int) * NM);
  }
}

void MatchFile::UpdateHeaderAndRecordLoc() {
  // the file _header
  if (_header_changed) {
    lseek(_fid, 0, SEEK_SET);
    // VERIFY(_header.version == MATCH_FILE_LATEST_VERSION);
    write(_fid, &_header, sizeof(_header));
    _header_changed = 0;
  }

  if (_header.definition_size > _recordloc_unchanged) {
    // recordloc
    lseek(_fid, sizeof(_header) + _recordloc_unchanged, SEEK_SET);
    write(_fid, &_record_definition[_recordloc_unchanged],
           _header.definition_size - _recordloc_unchanged);
    _recordloc_unchanged = _header.definition_size;
  }
}

///////////////////////
int MatchFile::GetPMatch(const char* image_path, int FC, int& NM,
                         Points<int>& matches) {
  int index;
  char relative_path[MAX_PATH];
  GetRelativePath(image_path, relative_path);
  index = GetImageIndex(relative_path);
  if (index < 0) {
    NM = 0;
    matches.resize(0, 0);
    return 0;
  } else {
    RecordLoc* loc = _match_records[index];

    if (loc->fcount == FC) {
      lseek(_fid, loc->read_loc, SEEK_SET);
      read(_fid, &NM, sizeof(int));
      lseek(_fid, sizeof(MatchRecordV3), SEEK_CUR);

      if (NM > 0 && NM <= FC) {
        matches.resize(NM, 2);
        read(_fid, matches[0], 2 * NM * sizeof(int));
      } else {
        matches.resize(0, 0);
        NM = 0;
      }
      return NM;
    } else {
      NM = 0;
      matches.resize(0, 0);
      return 0;
    }
  }
}

int MatchFile::GetPMatchR(const char* image_path, int FC, int& NM,
                          Points<int>& matches) {
  int index;
  char relative_path[MAX_PATH];
  GetRelativePath(image_path, relative_path);
  index = GetImageIndex(relative_path);
  if (index < 0) {
    NM = 0;
    matches.resize(0, 0);
    return 0;
  } else {
    RecordLoc* loc = _match_records[index];

    if (loc->fcount == FC) {
      lseek(_fid, loc->read_loc, SEEK_SET);
      read(_fid, &NM, sizeof(int));
      lseek(_fid, sizeof(MatchRecordV3), SEEK_CUR);

      if (NM > 0 && NM <= FC) {
        matches.resize(NM, 2);
        read(_fid, matches[1], NM * sizeof(int));
        read(_fid, matches[0], NM * sizeof(int));
      } else {
        matches.resize(0, 0);
      }
      return NM;
    } else {
      NM = 0;
      matches.resize(0, 0);
      return 0;
    }
  }
}

int MatchFile::GetMatchCount(const char* image_path, int& NM, int& NF) {
  int index;
  char relative_path[MAX_PATH];
  GetRelativePath(image_path, relative_path);
  index = GetImageIndex(relative_path);
  if (index < 0) {
    NM = NF = 0;
    return 0;
  } else {
    int counts[3];
    RecordLoc* loc = _match_records[index];
    lseek(_fid, loc->read_loc, SEEK_SET);
    read(_fid, counts, sizeof(counts));
    NM = counts[0];
    NF = counts[2];
    return 1;
  }
}

int MatchFile::GetIMatch(const char* image_path, int FC, TwoViewGeometry& tvg,
                         Points<int>& inliers) {
  int index;
  char relative_path[MAX_PATH];
  GetRelativePath(image_path, relative_path);
  index = GetImageIndex(relative_path);
  if (index < 0) {
    tvg.ResetGeometry();
    inliers.resize(0, 0);
    return 0;
  } else {
    RecordLoc* loc = _match_records[index];
    if (loc->fcount != FC) {
      tvg.ResetGeometry();
      inliers.resize(0, 0);
      printf("# features changed: [%s][%d->%d]!!!\r\n", image_path, loc->fcount,
             FC);
      return 0;
    } else {
      int NM;
      MatchRecordV3 rec;
      lseek(_fid, loc->read_loc, SEEK_SET);
      read(_fid, &NM, sizeof(int));
      read(_fid, &rec, sizeof(rec));
      if (NM > loc->fcount || NM < 0) {
        tvg.ResetGeometry();
        inliers.resize(0, 0);
        printf("ERROR: incorrect matching count [%d, %d]\r\n", NM, loc->fcount);
        return 0;

      } else if (NM == 0 || rec.version != MatchFile::MATCH_RECORD_V3 ||
                 rec.tvg.NF > loc->fcount) {
        tvg.ResetGeometry();
        inliers.resize(0, 0);
        // printf("#--BAD matching record---###\r\n");
        return 0;
      } else {
        //////////////////////////
        tvg.SetGeometry(rec.tvg);

        //
        if (tvg.NF > 0) {
          lseek(_fid, NM * 2 * sizeof(int), SEEK_CUR);
          inliers.resize(tvg.NF, 2);
          read(_fid, inliers[0], 2 * tvg.NF * sizeof(int));
        } else {
          inliers.resize(0, 0);
        }
        return NM;
      }
    }
  }
}

int MatchFile::GetIMatchR(const char* image_path, int FC, TwoViewGeometry& tvg,
                          Points<int>& inliers) {
  int index;
  char relative_path[MAX_PATH];
  GetRelativePath(image_path, relative_path);
  index = GetImageIndex(relative_path);
  if (index < 0) {
    tvg.ResetGeometry();
    inliers.resize(0, 0);
    return 0;
  } else {
    RecordLoc* loc = _match_records[index];
    if (loc->fcount != FC) {
      tvg.ResetGeometry();
      inliers.resize(0, 0);
      printf("# feature changed: [%s][%d->%d]!!!\r\n", image_path, loc->fcount,
             FC);
      return 0;
    } else {
      int NM;
      MatchRecordV3 rec;
      lseek(_fid, loc->read_loc, SEEK_SET);
      read(_fid, &NM, sizeof(int));
      read(_fid, &rec, sizeof(rec));
      if (NM > loc->fcount || NM < 0) {
        tvg.ResetGeometry();
        inliers.resize(0, 0);
        printf("ERROR: incorrect matching count [%d, %d]\r\n", NM, loc->fcount);
        return 0;

      } else if (NM == 0 || rec.version != MatchFile::MATCH_RECORD_V3 ||
                 rec.tvg.NF > loc->fcount) {
        tvg.ResetGeometry();
        inliers.resize(0, 0);
        // printf("#--BAD matching record---###\r\n");
        return 0;
      } else {
        ///
        tvg.SetGeometryR(rec.tvg);
        //
        if (tvg.NF > 0) {
          lseek(_fid, NM * 2 * sizeof(int), SEEK_CUR);
          inliers.resize(tvg.NF, 2);
          read(_fid, inliers[1], tvg.NF * sizeof(int));
          read(_fid, inliers[0], tvg.NF * sizeof(int));
        } else {
          inliers.resize(0, 0);
        }
        return NM;
      }
    }
  }
}

}  // namespace visual_sfm
