About wbpz files
 The wbpz file is an archive file of Ansys Workbench, and it should be used in Ansys Workbench environment.
 These were created by Ansys Workbench 2021R2. Please note that files cannot be opened with releases before 2021R2.
 In order to execute the simulation, "usermatLib.dll" in the folder "user_files" created by expanding the wbpz file is required. Please use Ansys Workbench after makeing the following settings.
  1. Copy "usermatLib.dll" to an arbitrary directory (such as "C:\AnsCustom" in the example).
  2. Set the path environment variable on Windows as follows:
     ANS_USER_PATH=C:\AnsCustom
 The result data is not included in "ansys_sphdd.wbpz" due to the large file size. You need to execute it to see the result.
