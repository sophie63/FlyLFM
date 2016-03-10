dir1 = getDirectory("Choose Text Files Directory ");

dir2 = getDirectory("Choose Destination Directory ");

file=File.openDialog("Select a File"); 


//setBatchMode(true);

//for (i=0; i<150; i++) {
for (i=0; i<150; i++) {
	t=i+1;
	Textfile="t"+t+".txt";
	filestring=File.openAsString(dir1+Textfile); 
	rows=split(filestring, "\n"); 
//	for (j=0; j<1202; j++) {
	for (j=0; j<1202; j++) {
		write(rows[j]);
		write(j);
		write(i);
		write(Textfile);
		open(file);
		run("Make Substack...", "slices=1-11 frames="+(i+1));
		selectWindow("melodic_IC-1.nii");
		run("16-bit");
		run("Multiply...", "value="+rows[j]+" stack");
		saveAs("Tiff", dir2+"C"+IJ.pad((i+1),3)+"_711ss1_T"+IJ.pad((j+1),4)+".tif");
		close("melodic_IC-1.nii");
		close("C"+IJ.pad((i+1),3)+"_711ss1_T"+IJ.pad((j+1),4)+".tif");
		close("melodic_IC.nii");
	}
}

//setBatchMode(false);
