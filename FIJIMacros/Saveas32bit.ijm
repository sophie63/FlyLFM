input = "/home/sophie/Desktop/713ss1partodor/";
output = "/home/sophie/Desktop/713partodor32/";
setBatchMode(true); 
list = getFileList(input);
for (i = 0; i < list.length; i++)
        action(input, output, list[i]);
setBatchMode(false);

function action(input, output, filename) {	
        open(input + filename);      
        run("32-bit");
        saveAs("Tiff", output + filename);
        close();
}
