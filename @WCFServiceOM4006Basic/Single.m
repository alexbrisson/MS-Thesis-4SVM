function Single(obj)
%Single(obj)
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
    'Single', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://tempuri.org/IWCFServiceOM4006Basic/Single', ...
    soapMessage);
parseSoapResponse(response);
