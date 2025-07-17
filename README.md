README: How to Set Up GMSH on Windows for Command Line Usage

============================================================
INTRODUCTION
============================================================

In order to run the library, Gmsh must first be installed
and configured correctly. Follow the steps below to install
Gmsh and enable command-line access.

============================================================
1. DOWNLOAD AND EXTRACT GMSH
============================================================

1. Go to the official Gmsh website:
   https://gmsh.info

2. Download the Windows version (e.g., gmsh-4.X.X-Windows64.zip).

3. Extract the .zip file to a folder, such as:
   C:\Users\YourName\Downloads\gmsh-4.X.X-Windows64

============================================================
2. ADD GMSH TO SYSTEM PATH
============================================================

1. Press Win + S and search for:
   "Edit the system environment variables"

2. In the System Properties window, click:
   "Environment Variables..."

3. Under "System variables", find and select:
   Path → then click "Edit"

4. Click "New" and add the full path to the folder
   containing gmsh.exe, e.g.:
   C:\Users\YourName\Downloads\gmsh-4.X.X-Windows64

5. Click OK to close all windows and apply changes.

============================================================
3. TEST GMSH IN COMMAND PROMPT
============================================================

1. Open a new Command Prompt window:
   Win + R → type: cmd → press Enter

2. Type:
   gmsh

3. The Gmsh GUI should open if everything is set up properly.

============================================================
4. OPTIONAL: RUN A GEO FILE FROM CMD
============================================================

If you have a .geo file and want to generate a 3D mesh:

   gmsh yourfile.geo -3

This will output mesh files in the same directory.

============================================================
NOTES
============================================================

- Make sure to restart Command Prompt after updating the PATH.
- You can use Gmsh both with the GUI and from scripts or terminal.
- For Python scripting with Gmsh, install the API via:
  
      pip install gmsh

