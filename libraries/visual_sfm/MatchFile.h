////////////////////////////////////////////////////////////////////////////
//	File:		MatchFile.h
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

#ifndef MATCH_FILE_H_INCLUDED
#define MATCH_FILE_H_INCLUDED
#pragma once

#include <cstring>
#include <vector>
#include <string>
using std::vector;
using std::string;
#include "points.h"

#define DISABLE_NAMED_MUTEX

namespace visual_sfm {

////////////////////////////////////////////////////////////////////////
//+++To read the number of matches of iamge1 and
//image2------------------------------
// Ues member funcation
// -----MatchFile::GetMatchCount--------------------------------
// int NM;		//number of putative matches
// int NF;		//number of inlier matches
// MatchFile mat(image1_full_path); //image1
// if(mat.IsValid() && mat.GetMatchCount(image2_full_path, NM, NF)) {/*we have
// it*/}

////////////////////////////////////////////////////////////////////////
//+++To read back putative
//matches----------------------------------------------------
// use member function
// ----MatchFile::GetPMatch----------------------------------------
// int FC;		//feature count of image2, input
// int NM;		//number of putative matches, output
// Points<int>	matches;// feature index of putative matches, output
// MatchFile mat(image1_full_path); //image1's match database
// if(mat.IsValid() && mat.GetPMatch(image2_full_path, FC, NM, matches)) {/* we
// have it*/};

//////////////////////////////////////////////////////////////////////////
//+++To write puattive
//matching----------------------------------------------------------
// Use static function
// ----MatchFile::WritePMatch---------------------------------------
// You can call MatchFile::WritePMatch(
//		char* image1_full_path_no_extension,
//		char* image2_full_path_no_extension,
//		int feature_count_of_image1,
//		int feature_count_of_image2,
//	    int number_of_matches,
//		Points<int>& matches)
//
//		You should call matches.resize(number_of_matches, 2) to form the
//matrix
//		matches[0][] are the indices of image1
//		matches[1][] are the indices of image2

#define MAX_PATH 512

struct TwoViewGeometry {
  int NF, NE;     // ni fundamental/essential matrix
  int NH, NH2;    // ni homogrphy
  float F[3][3];  // fundamental matrix
  float R[3][3];  // rotation
  float T[3];     // translation
  float F1, F2;   //
  float H[3][3];  //
  float GE;       // geometric error
  float AA;       // triangulation angle
  TwoViewGeometry() {
    NF = NE = NH = NH2 = 0;
    AA = F1 = F2 = 0.0;
    GE = 100.0f;
  }
  // simple geometry
  TwoViewGeometry(int Num, float FMat[3][3]) {
    NF = Num;
    NE = NH = NH2 = 0;
    AA = F1 = F2 = 0.0;
    GE = 100.0f;
    memcpy(F, FMat, sizeof(F));
  }
  void ResetGeometry() {
    NF = NE = NH = NH2 = 0;
    AA = F1 = F2 = 0.0;
    GE = 100.0f;
  }
  void ExchangeView() {
    float RT[3][3] = {{R[0][0], R[1][0], R[2][0]},
                      {R[0][1], R[1][1], R[2][1]},
                      {R[0][2], R[1][2], R[2][2]}};
    float TT[3] = {-R[0][0] * T[0] - R[1][0] * T[1] - R[2][0] * T[2],
                   -R[0][1] * T[0] - R[1][1] * T[1] - R[2][1] * T[2],
                   -R[0][2] * T[0] - R[1][2] * T[1] - R[2][2] * T[2]};
    float FT[3][3] = {{F[0][0], F[1][0], F[2][0]},
                      {F[0][1], F[1][1], F[2][1]},
                      {F[0][2], F[1][2], F[2][2]}};
    memcpy(F, FT, sizeof(F));
    memcpy(R, RT, sizeof(R));
    memcpy(T, TT, sizeof(T));
    float TF = F1;
    F1 = F2;
    F2 = TF;
  }
  static void TransposeMatrix33(const float M[3][3], float MT[3][3]) {
    MT[0][0] = M[0][0];
    MT[0][1] = M[1][0];
    MT[0][2] = M[2][0];
    MT[1][0] = M[0][1];
    MT[1][1] = M[1][1];
    MT[1][2] = M[2][1];
    MT[2][0] = M[0][2];
    MT[2][1] = M[1][2];
    MT[2][2] = M[2][2];
  }
  void SetGeometryR(const TwoViewGeometry& g) {
    NF = g.NF;
    NE = g.NE;
    NH = g.NH;
    NH2 = g.NH2;
    GE = g.GE;
    AA = g.AA;
    F1 = g.F2;
    F2 = g.F1;
    if (g.NF > 0) {
      TransposeMatrix33(g.F, F);
      if (g.NE > 0) {
        /////////////////////////////////////////////////////////////
        TransposeMatrix33(g.R, R);
        //////////////////////////////////////////////////////////////////
        T[0] = -R[0][0] * g.T[0] - R[0][1] * g.T[1] - R[0][2] * g.T[2];
        T[1] = -R[1][0] * g.T[0] - R[1][1] * g.T[1] - R[1][2] * g.T[2];
        T[2] = -R[2][0] * g.T[0] - R[2][1] * g.T[1] - R[2][2] * g.T[2];
      }
    }
    if (g.NH2 > 0) TransposeMatrix33(g.H, H);
  }

  void SetGeometry(const TwoViewGeometry& g) {
    NF = g.NF;
    NE = g.NE;
    NH = g.NH;
    NH2 = g.NH2;
    GE = g.GE;
    AA = g.AA;
    F1 = g.F1;
    F2 = g.F2;

    if (g.NF > 0) {
      memcpy(F, g.F, sizeof(F));
      if (g.NE > 0) {
        memcpy(R, g.R, sizeof(R));
        T[0] = g.T[0];
        T[1] = g.T[1];
        T[2] = g.T[2];
      }
    }
    if (g.NH2 > 0) memcpy(H, g.H, sizeof(H));
  }
};

class NamedMutex;

class MatchFile {
  enum {
    MATCH_FILE_VERSION_3 = ('M' + ('T' << 8) + ('0' << 16) + ('3' << 24)),
    MATCH_FILE_LATEST_VERSION = MATCH_FILE_VERSION_3,
    MATCH_RECORD_V3 = ('M' + ('R' << 8) + ('V' << 16) + ('3' << 24)),
  };

  struct MatchRecordV3 {
    int version;
    TwoViewGeometry tvg;
    int extra[7];  // reserved data for future change;
    ////////////////////////////////////////////////////
    MatchRecordV3() {
      version = MatchFile::MATCH_RECORD_V3;
      memset(extra, 0, sizeof(extra));
    }
  };

  // image file can change..
  // number of features can change.
  // exif focal length can change.
  struct FileHeader {
    int version;
    int file_count;
    int definition_size;  //
    int definition_buf;   //
    int feature_count;  // this might change..
  };

  struct RecordLoc {
    int fcount;
    int read_loc;    // read location
    int block_size;  // block size
    int trash_size;  // trash size after
    int extra_size;  //
    char file_name[4];
  };

 private:
  static int multi_thread_mode;
  static int file_title_mode;
  static int delay_header_update;
  static int record_reservation;
  static int num_match_verification;

 private:
  NamedMutex* _mutex;
  int _fid;
  int _opened_forwrite;
  int _header_changed;
  size_t _recordloc_unchanged;
  int _matchfile_mode;
  FileHeader _header;
  char _filepath[MAX_PATH];
  const char* _title;
  vector<char> _record_definition;
  vector<RecordLoc*> _match_records;

 private:
  void DeleteMatchRecord(const char* pattern = "extra");
  void GetRecordList();
  void MoveRecordToEnd(int i);
  void GetMatchFilePath(char* match_path);
  void GetMatchFolderPath(char* match_path);
  int AddImageMatch(const char* relative_path);
  int GetImageIndex(const char* relative_path);
  void GetRelativePath(const char* image_match, char* relative_path);
  void UpdateHeaderAndRecordLoc();
  void ResetMatchFile();

 public:
  // print out messages
  static int (*matchfile_printf)(const char* format, ...);
  static void SetPrintFunction(int (*printf_func)(const char* format, ...));
  static void GetRelativePath(const char* path1, const char* path2,
                              char* relative_path);
  static void SetMultiThreadMode(int multi_thread) {
    multi_thread_mode = multi_thread;
  }
  static void SetFileTitleMode(int use_title) { file_title_mode = use_title; }
  static void SetRecordReservation(int num) {
    record_reservation = (num > 1 ? num : 1);
  }

 public:
  MatchFile();
  MatchFile(const char* image_path, int guided = 0);
  ~MatchFile();
  int IsValid() { return _fid > 0; }
  int IsMatchFileOfImage(const char* fullpath, const char* title);

 public:
  int OpenMatchFile(const char* image_path, int write, int mode = 0);
  int MakeWritable();
  int HaveMatchRecord(const char* absolute_path);
  void DeleteMatchFile();
  void VerifyFeatureCount(int feature_count);
  void CloseMatchFile();
  void GetMatchedImageList(vector<string>& paths);
  void SaveSubsetMatch(vector<string>& paths, vector<int>& fc);

 public:
  //////////////////////////read matches
  int GetMatchCount(const char* image_path, int& NM, int& NF);
  int GetPMatch(const char* image_path, int FC, int& NM, Points<int>& matches);
  int GetIMatch(const char* image_path, int FC, TwoViewGeometry& tvg,
                Points<int>& inliers);
  int GetPMatchR(const char* image_path, int FC, int& NM, Points<int>& matches);
  int GetIMatchR(const char* image_path, int FC, TwoViewGeometry& tvg,
                 Points<int>& inliers);

 public:
  //////////////////////////////write matches
  void WriteIMatch(const char* image_match, int reverse, int NM,
                   TwoViewGeometry& tvg, Points<int>& inliers);
  void WritePMatch(const char* image_match, int FC, int NM,
                   Points<int>& matches, int reverse);

 public:
  ////////////////////static functions
  static void WriteIMatch(const char* image1, const char* image2, int NM,
                          TwoViewGeometry& tvg, Points<int>& inliers);
  static void WriteIMatch(MatchFile* mat, const char* image1,
                          const char* image2, int NM, TwoViewGeometry& tvg,
                          Points<int>& inliers);
  static void WritePMatch(const char* image1, const char* image2, int fc1,
                          int fc2, int NM, Points<int>& matches);
  static void WriteFMatch(const char* image1, const char* image2, int NF,
                          int* index1, int* index, float F[3][3]);
  ///////////////////////////////////////
  friend struct MatchFile::MatchRecordV3;
};

}  // namespace visual_sfm

#endif
