#include <gtest/gtest.h>
#include <rocksdb/db.h>

#include "theia/matching/image_pair_match.h"
#include "theia/matching/keypoints_and_descriptors.h"
#include "theia/matching/rocksdb_features_and_matches_database.h"

namespace theia {
namespace {
static std::string db_directory = THEIA_DATA_DIR + std::string("/database");

std::string RandomString(size_t length) {
  auto randchar = []() -> char {
    const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[rand() % max_index];
  };
  std::string str(length, 0);
  std::generate_n(str.begin(), length, randchar);
  return str;
}
}  // namespace

TEST(RocksDbFeaturesAndMatchesDatabase, Setup) {
  RocksDbFeaturesAndMatchesDatabase db(db_directory);
  rocksdb::DestroyDB(db_directory, rocksdb::Options());
}

TEST(RocksDbFeaturesAndMatchesDatabase, PutFeature) {
  static const std::string kImageName = "image_name";
  static const int kNumFeatures = 1000;

  // Create some features.
  KeypointsAndDescriptors features;
  features.keypoints.resize(kNumFeatures);
  features.descriptors.resize(kNumFeatures);
  for (int i = 0; i < kNumFeatures; i++) {
    features.keypoints[i] = Keypoint(i, i + 1, Keypoint::OTHER);
    features.descriptors[i].setRandom();
  }

  RocksDbFeaturesAndMatchesDatabase db(db_directory);

  // Add the features.
  db.PutFeatures(kImageName, features);

  // Get the features and ensure they are correct.
  const KeypointsAndDescriptors db_features = db.GetFeatures(kImageName);
  ASSERT_EQ(db_features.keypoints.size(), kNumFeatures);
  ASSERT_EQ(db_features.descriptors.size(), kNumFeatures);
  for (int i = 0; i < kNumFeatures; i++) {
    EXPECT_EQ(db_features.keypoints[i].x(), features.keypoints[i].x());
    EXPECT_EQ(db_features.keypoints[i].y(), features.keypoints[i].y());
    EXPECT_EQ(db_features.descriptors[i], features.descriptors[i]);
  }

  rocksdb::DestroyDB(db_directory, rocksdb::Options());
}

TEST(RocksDbFeaturesAndMatchesDatabase, GetFeatureFromInputDB) {
  static const std::string kImageName = "image_name";
  static const int kNumFeatures = 1000;

  // Create some features.
  KeypointsAndDescriptors features;
  features.keypoints.resize(kNumFeatures);
  features.descriptors.resize(kNumFeatures);
  for (int i = 0; i < kNumFeatures; i++) {
    features.keypoints[i] = Keypoint(i, i + 1, Keypoint::OTHER);
    features.descriptors[i].setRandom();
  }

  {
    RocksDbFeaturesAndMatchesDatabase db(db_directory);
    // Add the features.
    db.PutFeatures(kImageName, features);
    // Database closes when it goes out of scope.
  }

  {
    // Open the DB again and ensure it can retreive the features.
    RocksDbFeaturesAndMatchesDatabase db(db_directory);
    LOG(INFO) << "Num image features: " << db.NumImages();

    std::vector<std::string> image_names = db.ImageNamesOfFeatures();
    CHECK_GT(image_names.size(), 0);
    for (int i = 0; i < image_names.size(); i++) {
      LOG(INFO) << "Image name: " << image_names[i];
    }

    // Get the features and ensure they are correct.
    const KeypointsAndDescriptors db_features = db.GetFeatures(kImageName);
    ASSERT_EQ(db_features.keypoints.size(), kNumFeatures);
    ASSERT_EQ(db_features.descriptors.size(), kNumFeatures);
    for (int i = 0; i < kNumFeatures; i++) {
      EXPECT_EQ(db_features.keypoints[i].x(), features.keypoints[i].x());
      EXPECT_EQ(db_features.keypoints[i].y(), features.keypoints[i].y());
      EXPECT_EQ(db_features.descriptors[i], features.descriptors[i]);
    }
  }

  rocksdb::DestroyDB(db_directory, rocksdb::Options());
}

TEST(RocksDbFeaturesAndMatchesDatabase, ContainsFeature) {
  static const int kNumFeatures = 1000;
  static const int kStringLength = 64;

  // Open the DB again and ensure it can retreive the features.
  RocksDbFeaturesAndMatchesDatabase db(db_directory);

  std::vector<std::string> random_strings;
  for (int i = 0; i < kNumFeatures; i++) {
    random_strings.emplace_back(RandomString(kStringLength));
    db.PutFeatures(random_strings.back(), KeypointsAndDescriptors());
  }

  for (int i = 0; i < kNumFeatures; i++) {
    EXPECT_TRUE(db.ContainsFeatures(random_strings[i]));
  }

  rocksdb::DestroyDB(db_directory, rocksdb::Options());
}

TEST(RocksDbFeaturesAndMatchesDatabase, ImageNames) {
  static const int kNumFeatures = 1000;
  static const int kStringLength = 64;

  // Open the DB again and ensure it can retreive the features.
  RocksDbFeaturesAndMatchesDatabase db(db_directory);

  std::vector<std::string> random_strings;
  for (int i = 0; i < kNumFeatures; i++) {
    random_strings.emplace_back(RandomString(kStringLength));
    db.PutFeatures(random_strings.back(), KeypointsAndDescriptors());
  }
  std::sort(random_strings.begin(), random_strings.end());

  std::vector<std::string> image_names = db.ImageNamesOfFeatures();
  std::sort(image_names.begin(), image_names.end());

  ASSERT_EQ(random_strings.size(), image_names.size());
  for (int i = 0; i < random_strings.size(); i++) {
    EXPECT_EQ(random_strings[i], image_names[i]);
  }

  rocksdb::DestroyDB(db_directory, rocksdb::Options());
}

TEST(RocksDbFeaturesAndMatchesDatabase, PutMatch) {}

TEST(RocksDbFeaturesAndMatchesDatabase, GetMatchFromInputDB) {}

TEST(RocksDbFeaturesAndMatchesDatabase, ContainsMatch) {}
TEST(RocksDbFeaturesAndMatchesDatabase, MatchNames) {
  static const int kNumMatches = 1000;
  static const int kStringLength = 64;

  // Open the DB again and ensure it can retreive the features.
  RocksDbFeaturesAndMatchesDatabase db(db_directory);

  std::vector<std::pair<std::string, std::string>> random_strings;
  for (int i = 0; i < kNumMatches; i++) {
    random_strings.emplace_back(RandomString(kStringLength),
                                RandomString(kStringLength));
    db.PutImagePairMatch(
        random_strings[i].first, random_strings[i].second, ImagePairMatch());
  }
  std::sort(random_strings.begin(), random_strings.end());

  std::vector<std::pair<std::string, std::string>> match_names =
      db.ImageNamesOfMatches();
  std::sort(match_names.begin(), match_names.end());

  ASSERT_EQ(random_strings.size(), match_names.size());
  for (int i = 0; i < random_strings.size(); i++) {
    EXPECT_EQ(random_strings[i], match_names[i]);
  }

  rocksdb::DestroyDB(db_directory, rocksdb::Options());
}
}  // namespace theia
