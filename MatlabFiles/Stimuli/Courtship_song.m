%This script starts the excitation light with a LabJack and then the courstship song after
%5 sec (need a courtship song file)

clc %Clear the MATLAB command window
clear %Clear the MATLAB variables

ljasm = NET.addAssembly('LJUDDotNet'); %Make the UD .NET assembly visible in MATLAB
ljudObj = LabJack.LabJackUD.LJUD;
%Used for casting a value to a CHANNEL enum
chanType = LabJack.LabJackUD.CHANNEL.LOCALID.GetType;
[ljerror, ljhandle] = ljudObj.OpenLabJack(LabJack.LabJackUD.DEVICE.U3, LabJack.LabJackUD.CONNECTION.USB, '0', true, 0);

%illumination on
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 1, 0, 0);
ljudObj.GoOne(ljhandle);

pause(5)
load('C:\Users\Sophie Aimon\Desktop\droso_song_short.mat');
gong = audioplayer(data2*15, fs);
play(gong);

ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 0, 0, 0);

