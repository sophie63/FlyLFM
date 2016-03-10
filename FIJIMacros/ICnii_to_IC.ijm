
file=File.openDialog("Select a File"); 
dir2 = getDirectory("Choose Destination Directory ");


setBatchMode(true);

//for (i=0; i<150; i++) {
for (i=0; i<106; i++) {
		open(file);
		run("Make Substack...", "slices=1-32 frames="+(i+1));
		// selectWindow("melodic_IC-1.nii");
		run("16-bit");
		saveAs("Tiff", dir2+"C"+IJ.pad((i+1),3)+"_711ss1.tif");
		close("C"+IJ.pad((i+1),3)+"_711ss1.tif");
		close("melodic_IC.nii");
}

setBatchMode(false);
