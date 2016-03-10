
format = getFormat();
dir2 = getDirectory("Choose Destination Directory ");


//setBatchMode(true);

for (j=1;j<8;j++) {
	dir1 = "/home/sophie/Desktop/721ss2"+j+"/";

	list = getFileList(dir1);

	
 open(dir1+list[0]);
selectWindow(list[0]);
 open(dir1+list[1]);
selectWindow(list[1]);
print(dir1+list[0]);

run("Concatenate...", "  title=[Concatenated Stacks.tif] image1="+list[0]+" image2="+list[1]+" image3=[-- None --]");

close(list[0]);
close(list[1]);

for (i=2; i<list.length; i++) {
 showProgress(i+1, list.length);
 open(dir1+list[i]);
 selectWindow(list[i]);

 selectWindow("Concatenated Stacks.tif");
 run("Concatenate...", "  title=[Concatenated Stacks.tif] image1=[Concatenated Stacks.tif] image2="+list[i]+" image3=[-- None --]");

 close(list[i]);
}

saveAs(format, dir2+"721"+IJ.pad(j,1)+".tif");

close();

}

setBatchMode(false);
 
function getFormat() {
 formats = newArray("TIFF", "8-bit TIFF", "JPEG", "GIF", "PNG",
 "PGM", "BMP", "FITS", "Text Image", "ZIP", "Raw");
 Dialog.create("Batch Convert");
 Dialog.addChoice("Convert to: ", formats, "TIFF");
 Dialog.show();
 return Dialog.getChoice();
}




