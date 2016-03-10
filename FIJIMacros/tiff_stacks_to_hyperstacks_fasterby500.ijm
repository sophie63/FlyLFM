dir1 = getDirectory("Choose Source Directory ");
format = getFormat();
dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);


setBatchMode(true);

 open(dir1+list[0]);
selectWindow(list[0]);
 open(dir1+list[1]);
selectWindow(list[1]);

run("Concatenate...", "  title=[Concatenated Stacks.tif] image1="+list[0]+" image2="+list[1]+" image3=[-- None --]");

close(list[0])
close(list[1])

j=1;
for (i=2; i<list.length; i++) {
 showProgress(i+1, list.length);
 open(dir1+list[i]);
 selectWindow(list[i]);

 selectWindow("Concatenated Stacks.tif");
 run("Concatenate...", "  title=[Concatenated Stacks.tif] image1=[Concatenated Stacks.tif] image2="+list[i]+" image3=[-- None --]");

 close(list[i]);
print((i+1)%500);
print((i+1));

if ((i+1)%500 < 1) {
saveAs(format, dir2+j+".tif");
selectWindow(j+".tif");
close();
j=j+1;
print(j);
run("Concatenate...", "  title=[Concatenated Stacks.tif] image1="+list[0]+" image2="+list[1]+" image3=[-- None --]");
}

setBatchMode(false);

}
 
function getFormat() {
 formats = newArray("TIFF", "8-bit TIFF", "JPEG", "GIF", "PNG",
 "PGM", "BMP", "FITS", "Text Image", "ZIP", "Raw");
 Dialog.create("Batch Convert");
 Dialog.addChoice("Convert to: ", formats, "TIFF");
 Dialog.show();
 return Dialog.getChoice();
}




