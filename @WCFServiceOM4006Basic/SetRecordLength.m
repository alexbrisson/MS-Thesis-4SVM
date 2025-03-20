function SetRecordLength(obj,recordLength)
%SetRecordLength(obj,recordLength)
%
%     Input:
%       recordLength = (unsignedInt)
%   
%     Output:

% Build up the argument lists.
values = { ...
   recordLength, ...
   };
names = { ...
   'recordLength', ...
   };
types = { ...
   '{http://www.w3.org/2001/XMLSchema}unsignedInt', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://tempuri.org/', ...
    'SetRecordLength', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://tempuri.org/IWCFServiceOM4006Basic/SetRecordLength', ...
    soapMessage);
parseSoapResponse(response);
