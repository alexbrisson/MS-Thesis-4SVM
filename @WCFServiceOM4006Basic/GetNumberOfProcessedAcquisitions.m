function GetNumberOfProcessedAcquisitionsResult = GetNumberOfProcessedAcquisitions(obj)
%GetNumberOfProcessedAcquisitions(obj)
%
%     Input:
%   
%     Output:
%       GetNumberOfProcessedAcquisitionsResult = (unsignedInt)

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
    'GetNumberOfProcessedAcquisitions', ...
    values,names,types,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://tempuri.org/IWCFServiceOM4006Basic/GetNumberOfProcessedAcquisitions', ...
    soapMessage);
GetNumberOfProcessedAcquisitionsResult = parseSoapResponse(response);
