input = "/home/sophie/Desktop/863ss2/";
output = "/home/sophie/Desktop/863ss2c/";
setBatchMode(true); 
list = getFileList(input);
for (i = 0; i < list.length; i++)
        action(input, output, list[i]);
setBatchMode(false);

function action(input, output, filename) {	
        open(input + filename);     
	makeRectangle(14, 50, 183, 92);
        run("Crop");
	//run("Make Substack...", "  slices=10-42");
        saveAs("Tiff", output+"ss1"+IJ.pad((i+1),4)+".tif");
        close();
}
