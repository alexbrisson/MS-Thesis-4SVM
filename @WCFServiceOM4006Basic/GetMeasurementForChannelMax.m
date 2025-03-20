function GetMeasurementForChannelMaxResult = GetMeasurementForChannelMax(obj,channel,name)
%GetMeasurementForChannelMax(obj,channel,name)
%
%     Input:
%       channel = (string)
%       name = (string)
%   
%     Output:
%       GetMeasurementForChannelMaxResult = (double)

% Build up the argument lists.
values = { ...
   channel, ...
   name, ...
   };
names = { ...
   'channel', ...
   'name', ...
   };
types = { ...
   '{http://www.w3.org/2001/XMLSchema}string', ...
   '{http://www.w3.org/2001/XMLSchema}string', ...
   };

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://tempuri.org/', ...
    'GetMeasurementForChannelMax', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://tempuri.org/IWCFServiceOM4006Basic/GetMeasurementForChannelMax', ...
    soapMessage);
GetMeasurementForChannelMaxResult = parseSoapResponse(response);
