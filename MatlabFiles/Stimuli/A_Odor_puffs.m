% This script controls a valve with a LabJack to blow an odor onto the fly (3 pulses,)

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

 pause(10)
       
for (i=1:3)
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 4, 1, 0, 0);

%ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 0, 0, 0);
ljudObj.GoOne(ljhandle);

 c = cputime

pause(3)
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 4, 0, 0, 0);
%ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 1, 0, 0);
ljudObj.GoOne(ljhandle);
 c = cputime

pause(5)
end

ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 0, 0, 0);
