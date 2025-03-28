function SetBlockSize(obj,blockSize)
%SetBlockSize(obj,blockSize)
%
%     Input:
%       blockSize = (unsignedInt)
%   
%     Output:

% Build up the argument lists.
values = { ...
   blockSize, ...
   };
names = { ...
   'blockSize', ...
   };
types = { ...
   '{http://www.w3.org/2001/XMLSchema}unsignedInt', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://tempuri.org/', ...
    'SetBlockSize', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://tempuri.org/IWCFServiceOM4006Basic/SetBlockSize', ...
    soapMessage);
parseSoapResponse(response);
