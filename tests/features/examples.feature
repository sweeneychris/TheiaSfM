
Feature: we have example projects

Scenario: Run the build reconstruction demo application
  Given build executable build/bin/build_reconstruction exists
  And the file applications/build_reconstruction_flags.txt exists
  When I run the executable with arguments --flagfile=applications/build_reconstruction_flags.txt
  Then the command is successful


Scenario: Run the extract features demo application
  Given build executable build/bin/extract_features exists
  And the file data/image/img1.png exists
  When I run the executable with arguments --input_images=data/image/img1.png
  Then the command is successful
  And the file img1.png.features is created
