% This script gives periodic flashes of light (after 10 s, 3sec on, 3sec on, controlled by a LabJack)

clc %Clear the MATLAB command window
clear %Clear the MATLAB variables

ljasm = NET.addAssembly('LJUDDotNet'); %Make the UD .NET assembly visible in MATLAB
ljudObj = LabJack.LabJackUD.LJUD;
%Used for casting a value to a CHANNEL enum
chanType = LabJack.LabJackUD.CHANNEL.LOCALID.GetType;
[ljerror, ljhandle] = ljudObj.OpenLabJack(LabJack.LabJackUD.DEVICE.U3, LabJack.LabJackUD.CONNECTION.USB, '0', true, 0);

%excitation light on
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 1, 0, 0);
ljudObj.GoOne(ljhandle);
       
pause(10)
      
for (i=1:5)
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 6, 1, 0, 0);

%ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 0, 0, 0);
ljudObj.GoOne(ljhandle);

 c = clock

pause(3)
ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 6, 0, 0, 0);
%ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 1, 0, 0);
ljudObj.GoOne(ljhandle);
 c = clock;

pause(3)
end

ljudObj.AddRequest(ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, 7, 0, 0, 0);
