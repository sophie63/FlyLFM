input = "/media/sophie/New Volume/Data100156/";
output = "/media/sophie/New Volume/Data100156c/";
File.makeDirectory(output);
setBatchMode(true); 
list = getFileList(input);

f=File.open('/home/sophie/testFileOrder.txt');



for (i = 0; i < list.length; i++)
        action(input, output, list[i]);
setBatchMode(false);

function action(input, output, filename) {	
        //open(input + filename);     
	//makeRectangle(437, 0, 1033, 256);
        //run("Crop");
        print(f,list[i]);
	//run("Make Substack...", "  slices=10-42");
        //saveAs("Tiff", output+"100693ss1c-"+IJ.pad((i+1),5)+".tif");
        //saveAs("Tiff", output+"Data100156_"+i+".tif");
        //close();
}
File.close('/home/sophie/testFileOrder.txt');