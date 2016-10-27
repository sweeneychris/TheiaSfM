.. _chapter-contributing:

=====================
Contributing to Theia
=====================

We welcome and encourage contributions to Theia, whether they are new features,
bug fixes or tests. The `Theia mailing list
<http://groups.google.com/group/theia-vision-library>`_ is the best place for
all development related discussions. Please consider joining it. If you have an
idea for how you'd like to contribute to Theia, please consider emailing the
list first to voice your idea. We can help you fine-tune your idea and this will
also help avoid duplicate work by somebody else who may be working on the same
feature.

Style and Testing
=================

We follow Google's `C++ Style Guide
<https://google.github.io/styleguide/cppguide.html>`_ and use git for version
control. We use `Gerrit <https://code.google.com/p/gerrit/>`_ to review changes
before commits. The git repository has been set up with
`GerritHub <http://gerrithub.io/>`_ so that the repository lives in GitHub, but
all reviews are performed with Gerrit.

When contributing substantial new code or bug fixes, please add unit tests to
ensure the usage of the code (or to prove the bug is fixed!).

CMake
=====

We use CMake to generate makefiles for Theia to maximize the cross-platform
usability. If you need to add a new library or a new file to Theia, you will
likely need to add that file to the CMakeLists.txt in src/theia (along with a
unit test!).

Developing for Theia
====================

Much of the instructions that follow in this section were borrowed and modified
from the `Ceres Solver
<http://homes.cs.washington.edu/~sagarwal/ceres-solver/stable/contributing.html>`_
project.

1. Download and configure ``git``.

   * Mac ``brew install git``.
   * Linux ``sudo apt-get install git``.
   * Windows. Download `msysgit
     <https://code.google.com/p/msysgit/>`_, which includes a minimal
     `Cygwin <http://www.cygwin.com/>`_ install.

2. Sign up for a `GitHub <http://github.com>`_ account. Accounts are free to register.

3. Clone the Theia ``git`` repository from GerritHub.

   .. code-block:: bash

      git clone https://review.gerrithub.io/sweeneychris/TheiaSfM


4. Build Theia, following the instructions in
   :ref:`chapter-building`.

   On Mac and Linux, the ``CMake`` build will download and enable
   the Gerrit pre-commit hook automatically. This pre-submit hook
   creates `Change-Id: ...` lines in your commits.

   If this does not work OR you are on Windows, execute the
   following in the root directory of the local ``git`` repository:

   .. code-block:: bash

      curl -o .git/hooks/commit-msg http://www.theia-sfm.org/_static/commit-msg
      chmod +x .git/hooks/commit-msg

5. Configure your GerritHub password.

  Sign into `https://review.gerrithub.io <https://review.gerrithub.io>`_, go to
  ``Settings`` then ``HTTP Password``. If no password exists, select ``Generate
  Password``. The username and password listed on this page will need to be
  entered when pushing your changes to the repo for review (see instructions
  below).

Submitting a change to Theia
============================

1. Make your changes against master or whatever branch you like. Commit your
   changes as one patch. This is critical for the review process; if you submit
   a change with multiple commits, each commit will become a separate review in
   Gerrit.

2. Push your changes to the Theia GerritHub instance:

   .. code-block:: bash

      git push origin HEAD:refs/for/master

   You will likely have to enter your GerritHub username and password. When the
   push succeeds, the console will display a URL showing the address of the
   review. Go to the URL and add reviewers; at this point this is only Chris.

3. Wait for a review.

4. Once review comments come in, address them. Please reply to each
   comment in Gerrit, which makes the re-review process easier. After
   modifying the code in your ``git`` instance, *don't make a new
   commit*. Instead, update the last commit using a command like the
   following:

   .. code-block:: bash

      git commit --amend -a

   This will update the last commit, so that it has both the original
   patch and your updates as a single commit. You will have a chance
   to edit the commit message as well. Push the new commit to Gerrit
   as before.

   Gerrit will use the ``Change-Id:`` to match the previous commit
   with the new one. The review interface retains your original patch,
   but also shows the new patch.

   Publish your responses to the comments, and wait for a new round
   of reviews.

5. Before submitting, make sure you are synced to the latest commit in the
   repo. To do this, simply run the command:

   .. code-block:: bash

      git pull --rebase origin master

   This will pull the latest changes without interfering with your current
   patch.
