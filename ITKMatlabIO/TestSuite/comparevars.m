%COMPAREVARS compare variables to within a given tolerance
%
%  Checks for variable comparison - returns true or false
%
%      flag = COMPAREVARS ( A, B )
%
%      flag = COMPAREVARS ( A, B, DCP, NaNCheck )
%              DCP      : number of decimal places to check (1e-8 default)
%              NaNCheck : if true (default) NaNs are treated as being equal
%
%                          entering [] in either DCP or NaNCheck will
%                          result in the code using the defaults.
%
%      flag = COMPAREVARS ( A, B, DCP, NaNCheck, debug )
%
%              debug    : flag to write check log to screen (default true)
%
% Example
%
%  %  create a structure with lots of difference data types
%   A.data.level.x = [1 2 3];
%   A.data.level(2).x = [1 2 3];
%   A.data.level(2).func = str2func ( 'plot' );
%   A.data.cell = cell(20,1);
%   A.other.log.flag = true;
%   A.dummyData.integers.test = int8(1);
%   A.dummyData.char(20).example = 'ABCDeFGH';
%   A.dummyData.img =i;
%
%   % change value x
%   B = A;
%   B.data.level(2).x = zeros(3,2,4,3);
%   flag = comparevars ( A, B )
%
%    see also isequal isequalwithequalnans test_comparevars
%

% Author    : Robert Cumming
% Copyright : Consultant Engineer Ltd
% email     : rcumming@matpi.com
% web       : www.matpi.com
%  $Id: comparevars.m 11 2014-04-08 18:17:54Z Bob $\n', mfilename );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function output = comparevars ( A, B, tolerance, nanCheck, debug )
  % validate which SVN the code is from
  if nargin == 1 && ischar ( A )                          % version information
    if strcmp ( A, '-version' )
      fprintf ( '%s Version: $Id: comparevars.m 11 2014-04-08 18:17:54Z Bob $\n', mfilename );
      return
    end
  end
  if nargin < 5; debug = true; end                                % debug - write out where the outputs are not equal
  if nargin < 4 || isempty ( nanCheck ); nanCheck = false; end                             % whether to
  if nargin < 3 || isempty ( tolerance )                            % check for defaults
    tolerance = 1e-8;                                              % set defaults
  end
  if iscell ( A )
    output = CheckAlmostEqualOnVarClass ( A, B, tolerance, nanCheck, debug, inputname(1), inputname(2) );    % if cell - go to check by type to recursively call
  elseif isstruct ( A )
    output = LoopAroundStructCheckEqual ( A, B, tolerance, nanCheck, debug, inputname(1), inputname(2) );    % is struct go to struct checking routine
  else
    output = CheckAlmostEqualOnVarClass ( A, B, tolerance, nanCheck, debug, inputname(1), inputname(2) );
  end
end
%---------------------------------------------------------------------------------------------------------------------
% Loop around structures to compare.
%---------------------------------------------------------------------------------------------------------------------
function output = LoopAroundStructCheckEqual ( A, B, tolerance, nanCheck, debug, path1, path2 )
  nA = length ( A );
  nB = length ( B );
  if isequal ( nA, nB )
    output = 1;
    for ii=1:length(A)
      subPath1 = sprintf ( '%s(%i)', path1, ii );  % build the subPath var which is passed into the child function
      subPath2 = sprintf ( '%s(%i)', path2, ii );  % build the subPath var which is passed into the child function
      % for each element (ii) of the structure call the struct function
      thisOutput = CheckAlmostEqualStruct ( A(ii), B(ii), tolerance, nanCheck, debug, subPath1, subPath2 );
      if output % if output true - copy the sub field output to it.  If output = false -> leave false.
        output = thisOutput;
      end
    end
  else
    if ( debug );
      fprintf ( 2, 'Check failed on length of structs - %s (%i) and %s (%i)\n', path1, nA, path2, nB );
    end
    output = false;
  end
end
%---------------------------------------------------------------------------------------------------------------------
% Method which checks structs are eqaul
%---------------------------------------------------------------------------------------------------------------------
function output = CheckAlmostEqualStruct ( A, B, tolerance, nanCheck, debug, path1, path2 )
  fnames1 = sort ( fieldnames ( A ));                                                 % Extract fieldnames from strucutre A
  fnames2 = sort ( fieldnames ( B ));                                                 % Extract fieldnmaes from structure B
  if isequal ( fnames1, fnames2 )                                                     % Check they are equal
    output = 1; % Init (in case length(fnames1) == 0 output would not be set)
    for ii=1:length(fnames1)                                                           % Loop for each sub field
      subPath1 = sprintf ( '%s.%s', path1, fnames1{ii} );
      subPath2 = sprintf ( '%s.%s', path2, fnames1{ii} );
      thisOutput = CheckAlmostEqualOnVarClass ( A.(fnames1{ii}), B.(fnames2{ii}), tolerance, nanCheck, debug, subPath1, subPath2 );  % Check on the type
      if output % if output true - copy the sub field output to it.  If output = false -> leave false.
        output = thisOutput;
      end
    end                                                                               %
  else                                                                                %
    output = false;                                                                   % If fields not equal - return false
    fnames1 = ExplodeCellToStr ( fnames1 );
    fnames2 = ExplodeCellToStr ( fnames2 );
    if ( debug );
      fprintf ( 2, 'Check failed on non equal fields in struct\n   %s contains - %s\n   %s contains - %s\n', path1, fnames1, path2, fnames2 );
    end
  end
end
%---------------------------------------------------------------------------------------------------------------------
% method to check all class types
%---------------------------------------------------------------------------------------------------------------------
function output = CheckAlmostEqualOnVarClass ( A, B, tolerance, nanCheck, debug, path1, path2 )
  % Check to see if they are both empty
  if isempty ( A ) && isempty ( B )
    output = 1;
  % Check for different class types
  elseif ~isequal ( class ( A ), class ( B ) )
    output = false;
    if ( debug );
      fprintf ( 2, 'Check failed different class types\n   %s (%s)\n   %s (%s)\n', path1, class ( A ), path2, class ( B ) );
    end
  else
    % Switch the type of variable
    switch class ( A );
      case 'cell'
        % check the length of the cells are the same
        s1 = size(A);
        s2 = size(B);
        if isequal ( s1, s2 )
          % loop for each cell
          output = 1;
          for ii=1:length(A)
            subPath1 = sprintf ( '%s{%i}', path1, ii );  % build the path string
            subPath2 = sprintf ( '%s{%i}', path2, ii );  % build the path string
            % if struct call the struct sub function - otherwise self call this function
            if isstruct (A{ii})
              thisOutput = CheckAlmostEqualStruct ( A{ii}, B{ii}, tolerance, nanCheck, debug, subPath1, subPath2 );
            else
              thisOutput = CheckAlmostEqualOnVarClass ( A{ii}, B{ii}, tolerance, nanCheck, debug, subPath1, subPath2 );
            end
            if output % if output true - copy the sub field output to it.  If output = false -> leave false.
              output = thisOutput;
            end
          end
        else
          % set output to false and return
          output = false;
          if ( debug );
            s1str = ConvertSizeVectorToString ( s1 );
            s2str = ConvertSizeVectorToString ( s2 );
            fprintf ( 2, 'Check failed on cell lengths \n   %s (%s)\n   %s (%s)\n', path1, s1str, path2, s2str );
          end
%           return
        end
      case 'char'  % check for char types - simple if equal - otherwise fail
        output = isequal ( A, B );
        if ( output == false && debug );
          fprintf ( 2, 'Check failed on char\n   %s = "%s"\n   %s = "%s"\n', path1, A, path2, B );
        end
      case 'struct'
        % call the loop for the elements of the structure  ii.e. var(1).X, var(2).X etc...
        s1 = length ( A );
        s2 = length ( B );
        if ~isequal ( s1, s2 )
          if ( debug );
            s1str = ConvertSizeVectorToString ( s1 );
            s2str = ConvertSizeVectorToString ( s2 );
            fprintf ( 2, 'Check failed on size of structs not the same\n   %s (%s)\n   %s (%s)\n', path1, s1str, path2, s2str );
          end
          output = false;
        else
          output = LoopAroundStructCheckEqual ( A, B, tolerance, nanCheck, debug, path1, path2 );
        end
      case { 'double' 'single' }
        % the main tolerance checking function - arrays, matrix etc....
        s1 = size(A);
        s2 = size(B);
        if ~isequal ( s1, s2 )
          s1str = ConvertSizeVectorToString ( s1 );
          s2str = ConvertSizeVectorToString ( s2 );
          if ( debug );
            fprintf ( 2, 'Check failed on size of variables not the same\n   %s (%s)\n   %s (%s)\n', path1, s1str, path2, s2str );
          end
          output = false;
        else
          A = A(:); % convert to column vector
          B = B(:); % convert to column vector
          % check to see if there are any NaN in the matrix
          % if there is remove the NaN from the array and self call to check the tolerance.
          if nanCheck && max ( isnan ( A ) ) || max ( isnan ( A ) )
            indexA = isnan ( A );
            indexB = isnan ( B );
            output = CheckAlmostEqualOnVarClass ( A(~indexA), B(~indexB), tolerance, nanCheck, debug, path );
          else
            % check for real/img differences between the variables.
            imgCheck(2) = isreal ( B );
            imgCheck(1) = isreal ( A );
            if mod ( sum ( imgCheck ), 2 )
              if ( debug );
                if imgCheck(1)
                  fprintf ( 2, 'Check failed on real/img comparison\n   %s (real)\n   %s (img)\n', path1, path2 );
                else
                  fprintf ( 2, 'Check failed on real/img comparison\n   %s (img)\n   %s (real)\n', path1, path2 );
                end
              end
              output = false;
            else
              % is tolerance is 0 - use inbuilt matlab check isequal
              if tolerance == 0
                if nanCheck  % should never get here, due to the if part of this - but put it in case.
                  output = isequaln ( A, B );
                else
                  output = isequal ( A, B );
                end
              else
                % Check the max of the absolute difference.
                output = max(abs(A-B))<tolerance;
              end
              if ( output == false && debug );
                % print to screen if requested.
                fprintf ( 2, 'Check failed on tolerance (%i) check\n   %s\n   %s\n', tolerance, path1, path2 );
              end
            end
          end
        end
      case 'logical'
        % check for logical values
        output = isequal ( A, B );                              % Check that logical are equal
        if ( output == false && debug );
          fprintf ( 2, 'Check failed on logical\n   %s=%i\n   %s=%i\n', path1, A, path2, B );
        end
      case 'function_handle'                                               %
        % check for function handles
        output = isequal ( A, B );                              % Check that function handles are equal
        if ( output == false && debug );
          fprintf ( 2, 'Check failed on function handle\n   %s="%s"\n   %s="%s"\n', path1, func2str(A), path2, func2str(B) );
        end
      case { 'int8' 'int16' 'int32' 'int64' }
        % check for ints
        A = A(:); % convert to array
        B = B(:); % convert to array
        output = isequal ( A, B );                              % Check that ints are equal
        output = max(abs(output));
        if ( output == false && debug );
          fprintf ( 2, 'Check failed on %s check - values different\n   %s\n   %s\n', class ( A ), path1, path2 );
        end
      otherwise
        % Otherwise write unsupported class
        %  and perform the normal isequal function
        fprintf ( 2, 'comparevars:Unsupported class %s\n', class(A) );
        output = isequal ( A, B );
        if ( output == false && debug );
          fprintf ( 2, 'Check failed on unsupported check\n   %s\n   %s\n', path1, path2 );
        end
    end
  end
end
% convert a vector into a formatted string.
function str = ConvertSizeVectorToString ( s1 )
  str = '';
  for ii=1:length(s1);
    str = sprintf ( '%s%i,', str, s1(ii) );
  end
  try; str(end) = []; end;
end
% convert a cell of char into a string.
function str = ExplodeCellToStr ( fnames )
  str = '';
  for ii=1:length(fnames)
     str = sprintf ( '%s%s,', str, fnames{ii} );
  end
  try; str(end) = []; end;
end
