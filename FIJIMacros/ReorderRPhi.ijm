input = "/home/sophie/Downloads/flow/";
output = "/home/sophie/Downloads/flow2/";
setBatchMode(true); 
list = getFileList(input);

for (i = 0; i < list.length; i++)
        action(input, output, list[i]);
        print(i)
setBatchMode(false);

function action(input, output, filename) {	
        open(input + filename); 
        print(input+filename);   
        run("Re-order Hyperstack ...", "channels=[Slices (z)] slices=[Channels (c)] frames=[Frames (t)]");
        run("NIfTI-1", "save="+output + filename);
        close();
}
