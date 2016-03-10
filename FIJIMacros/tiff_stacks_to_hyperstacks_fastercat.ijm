dir = getDirectory("Choose Source Directory ");
format = getFormat();

for (j=1; j<26; j++) {

dir1=dir+j+"/";
list = getFileList(dir1);


setBatchMode(true);

 open(dir1+list[0]);
selectWindow(list[0]);
 open(dir1+list[1]);
selectWindow(list[1]);

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

selectWindow("Concatenated Stacks.tif");

name=dir+j+".tif";
print(name);
saveAs(format, name);

close();
}

setBatchMode(false);

//for (j=1; j<9; j++) {
//name=dir+j+".tif";
//open(name);
//}

//run("Concatenate...", "  title=ss2.tif image1=1.tif image2=2.tif image3=3.tif image4=4.tif image5=5.tif image6=6.tif image7=7.tif image8=8.tif image9=[-- None --]");

 
function getFormat() {
 formats = newArray("TIFF", "8-bit TIFF", "JPEG", "GIF", "PNG",
 "PGM", "BMP", "FITS", "Text Image", "ZIP", "Raw");
 Dialog.create("Batch Convert");
 Dialog.addChoice("Convert to: ", formats, "TIFF");
 Dialog.show();
 return Dialog.getChoice();
}




