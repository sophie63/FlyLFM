
dir2 = getDirectory("Choose Destination Directory ");
file=File.openDialog("Select a File"); 

open(file);

for (i=1; i<=381; i++) {
	run("Make Substack...", "slices=1-34 frames="+i);
	saveAs("Tiff", dir2+"713ss1partodor"+IJ.pad(i,3)+".tif");
}

