function LoadPresets(obj,fileName)
%LoadPresets(obj,fileName)
%
%     Input:
%       fileName = (string)
%   
%     Output:

% Build up the argument lists.
values = { ...
   fileName, ...
   };
names = { ...
   'fileName', ...
   };
types = { ...
   '{http://www.w3.org/2001/XMLSchema}string', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://tempuri.org/', ...
    'LoadPresets', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://tempuri.org/IWCFServiceOM4006Basic/LoadPresets', ...
    soapMessage);
parseSoapResponse(response);
