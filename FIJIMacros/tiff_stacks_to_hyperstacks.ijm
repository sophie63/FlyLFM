dir1 = getDirectory("Choose Source Directory ");
format = getFormat();
dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);


setBatchMode(true);


 open(dir1+list[0]);
selectWindow(list[0]);
 open(dir1+list[1]);
selectWindow(list[1]);

run("Concatenate...", "  title=[Concatenated Stacks] image1="+list[0]+" image2="+list[1]+" image3=[-- None --]");
 saveAs(dir2+"Concatenated Stacks.tif");
close()

for (i=2; i<list.length; i++) {
 showProgress(i+1, list.length);
 open(dir1+list[i]);
 selectWindow(list[i]);
 open(dir2+"Concatenated Stacks.tif");
 selectWindow("Concatenated Stacks.tif");
 run("Concatenate...", "  title=[Concatenated Stacks] image1=[Concatenated Stacks.tif] image2="+list[i]+" image3=[-- None --]");
  saveAs(dir2+"Concatenated Stacks.tif");
 close();
}

open(dir2+"Concatenated Stacks.tif");
 selectWindow("Concatenated Stacks.tif");
run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=20 frames="+list.length+" display=Grayscale");

 saveAs(format, dir2+list[1]);

close()

setBatchMode(false);
 
function getFormat() {
 formats = newArray("TIFF", "8-bit TIFF", "JPEG", "GIF", "PNG",
 "PGM", "BMP", "FITS", "Text Image", "ZIP", "Raw");
 Dialog.create("Batch Convert");
 Dialog.addChoice("Convert to: ", formats, "TIFF");
 Dialog.show();
 return Dialog.getChoice();
}




