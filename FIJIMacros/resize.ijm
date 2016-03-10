input = "/home/sophie/Desktop/825ss2c/";
output = "/home/sophie/Desktop/825ss1c/";
setBatchMode(true); 
list = getFileList(input);
for (i = 0; i < list.length; i++)
        action(input, output, list[i]);
setBatchMode(false);

function action(input, output, filename) {	
        open(input + filename);   
        run("Resize ", "sizex=86.0 sizey=61.0 method=Least-Squares interpolation=Cubic unitpixelx=true unitpixely=true");
	run("Make Substack...", "  slices=1-33-2");
        saveAs("Tiff", output+"ss1"+IJ.pad((i+1),4)+".tif");
        close();
}
