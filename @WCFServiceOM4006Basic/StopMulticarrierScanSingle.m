function StopMulticarrierScanSingle(obj)
%StopMulticarrierScanSingle(obj)
%
%     Input:
%   
%     Output:

% Build up the argument lists.
values = { ...
   };
names = { ...
   };
types = { ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://tempuri.org/', ...
    'StopMulticarrierScanSingle', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://tempuri.org/IWCFServiceOM4006Basic/StopMulticarrierScanSingle', ...
    soapMessage);
parseSoapResponse(response);
