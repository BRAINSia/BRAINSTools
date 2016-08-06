classdef MatlabITKImageObject < handle
    % MATLABITKIMAGEobjHdlECT is a matlab image class with similar
    %   functionality as an ITK image objHdlect.
    %
    % This class is used to write/prototype algorithms in matlab using
    %   ITK structure so that the algorithm can be easily translated into
    %   c++/ITK and used for research.
    % references
    % webpages
    % handle makes this objHdlect pass by simulated reference
    % Write the presentation first

    %TODO - make these attributes private
    properties (Access = public)  %% ITK image attributes
        m_ImageDimension;%template parameter in ITK %TODO make this private only allowed to modify in the constructor
        m_Spacing;
        m_Origin;
        m_Direction;
        m_InverseDirection;
        m_IndexToPhysicalPoint;
        m_PhysicalPointToIndex;%n+1xn+1 matrix wheere n is the image dimension
        %TODO - remove m_Size, GetLargestPossibleRegionSize or GetRegionSize() note: this is convenience function
        % itk equivalent to im->getlargestpossible region().getsize()
        %its redundant with introspection of m_Buffer,
        %which is a 3D vector ( always verify m_Size is always == to
        %size(m_Buffer) )
        m_Size; %largest possible region size, does not support streaming
        %the matlab implemetnation only uses m_Size for defining the image,
        %iTK has a more expressive definition of a region that allows for
        %streaming and buffering. These region capabilities are NOT
        %available in this matlab interface
       % m_LargestPossibleRegion;% size and start index are in here
        %m_RequestedRegion;%unneeded
        %m_BufferedRegion;%unneeded explicitly not implemented for this interface
        m_Buffer;
    end

    %nrrd (DWI)attributes from nrrdCommon.h
    properties (Access = public)
        % documentation available at
        % http://teem.sourceforge.net/nrrd/lib.html
        % from nrrd matlab struct
        data;%/* the data in memory */
     %redundant image dimension above -REMOVE
        space;%/* from nrrdSpace* enum, and often implies the value of spaceDim */
     %REMOVE SPACE DIRECTION, orientation * spacing
        spacedirections;%/* the vector, in "space" (as described by
                        %nrrd->space and/or nrrd->spaceDim), from one
                        %sample to the next sample along this axis.  It
                        %is the column vector of the transform from
                        %index space to "space" space */
        centerings;%in ITK it is always centered in the middle of each cell
        kinds;%   ** For describing the information along one axis of an array.  This is
%                 ** most important for clarifying the representation of non-scalar
%                 ** data, in order to distinguish between axes that are genuine image
%                 ** domain axes, and axes that exist just to store the multiple
%                 ** attributes per sample.  One could argue that this information
%                 ** should be per-array and not per-axis, but you still have to
%                 ** indicate which one of the axes is the attribute axis.  And, if you
%                 ** have, say, the gradient of RGB colors, you want the per-pixel 3x3
%                 ** array to have those two attribute axes tagged accordingly.
        spaceunits;%/* units for coordinates of space */
    %if not LPS everthing else is tossed
    %validation WRARNING and check if data is not LPS
    %ITK cannot return anything other than LPS
    %
        spacedefinition;%** Identifies the space in which which the origin and direction
% **                      vectors have their coordinates measured.  When a direction is named
% **                      here (like "Left" or "Anterior"), that implies a basis vector that
% **                      points in that direction, along which that coordinate becomes *larger*
% **                      (this is the opposite of MetaIO, for example).
    %TODO Remove spaceorigin redundant
        spaceorigin;%/* the location of the center the first
                    %(lowest memory address) array sample,
                    %regardless of node-vs-cell centering */
        measurementframe;%/* if spaceDim is non-zero, this may store
%                            a spaceDim-by-spaceDim matrix which
%                            transforms vector/matrix coefficients
%                            in the "measurement frame" to those in
%                            the world space described by spaceDim
%                            (and hopefully space).  Coeff [i][j] is
%                             *column* i & *row* j, which is probably
%                            the *transpose* of what you expect.
%                            There are no semantics linking this to
%                            the "kind" of any axis, for a variety
%                            of reasons */
        modality;
        bvalue;
        gradientdirections;
    end

    methods (Access = public)

        % Constructor for default
        % Constructor for provided image dimension only
        % for Copy Construction
        % for construction from MATLAB NRRD Struct
        % for construction from NRRD filename/path
        % matlab OOP overloading
        function objHdl = MatlabITKImageObject(rhs)
            % default constructor
            if nargin == 0
                objHdl.SetImageDimension(3);
                DEBUG( 'Constructed from Default Constructor' )

            % default with image dimension
            elseif isa(rhs,'uint64') %TODO integer 64 bit only!!!!!!!!!!!
                objHdl.SetImageDimension(rhs);
                DEBUG( 'Constructed from Default Constructor with Specified Image Dimension' )

            % Copy Construction
            elseif  isa(rhs,'MatlabITKImageObject')
                fns = properties(rhs);
                for i=1:length(fns)
                    objHdl.(fns{i}) = rhs.(fns{i});
                end
                DEBUG( 'Constructed from Copy Constructor' )

            %construction from some nrrd file
            else
                % Construction from NRRD Struct
                if isstruct(rhs)
                    DEBUG( 'Constructed from NRRD Struct' )

                % Construction from NRRD Filename
                else
                    try
                        rhs = itkLoadWithMetadata(rhs);
                    catch ME
                        disp(ME)
                        disp( 'currently unsupported by itkLoadWithMetadata' )
                        try
                            disp( 'using nrrdLoadWithMetadata' )
                            rhs = nrrdLoadWithMetaData(rhs);
                        catch ME1
                            disp(ME1)
                            disp( 'ERROR - unsupported by itkLoad and nrrdLoad' )
                        end
                    end
                    DEBUG( 'Constructed from NRRD Filename' )
                end
                % get nrrd specific information and data and assign to
                % properties
                nrrdNames = fieldnames(rhs);
                for i = 1 : length(nrrdNames)
                    objHdl.(nrrdNames{i}) = rhs.(nrrdNames{i});
                end
                % assign ITK specific properties from nrrd metadata
                objHdl.SetImageDimension( rhs.space );
                objHdl.SetOrigin( rhs.spaceorigin );
                objHdl.SetDirection( rhs.spacedirections );
            end
        end%constructor

        % Allocate space in each ITK property for the specified
        % ImageDimension
        % TODO ONLY BE CALLED BY CONSTRUCTOR
        % ACCESS PRIVATE
        function SetImageDimension( objHdl, ImageDimension )
            objHdl.m_ImageDimension = ImageDimension;
            %TODO do this somewhere else
            objHdl.m_Spacing        = ones([ImageDimension,1]);
            objHdl.m_Origin         = ones([ImageDimension,1]);
            objHdl.m_Direction      = eye(ImageDimension);
            objHdl.m_InverseDirection       = eye(ImageDimension);
            objHdl.m_IndexToPhysicalPoint   = ones([ImageDimension,1]);
            objHdl.m_PhysicalPointToIndex   = ones([ImageDimension,1]);
        end

        % Print all properties of the Object
        function PrintSelf(objHdl)
            props = properties(objHdl);
            for i=1:length(props)
                % suppress data and gradientdirections output
                if ~strcmp(props(i),'data') &&...
                        ~strcmp(props(i),'gradientdirections')
                    disp( props{i} )
                    disp( objHdl.(props{i}) )
                end
            end
        end

        %------------SETTERS-----------------------------------------------
        % set the origin
        function SetOrigin(objHdl,origin)
            objHdl.m_Origin = origin;
        end

        % set the direction cosines matrix
        % calculate the inverse direction cosines matrix
        % calculate the index to physical point conversion matrix
        % calculate the physical point to index conversion matrix
        function SetDirection( objHdl, direction )
            if objHdl.m_Direction ~= direction
                objHdl.m_Direction = direction;
                objHdl.m_InverseDirection = direction^(-1);
                objHdl.m_IndexToPhysicalPoint = ...%index to physical point is 3x4 with origin values
                    objHdl.m_Direction * objHdl.m_Spacing;% + origin;%padd these 4x4 with zeros so spacing becomes
                objHdl.m_PhysicalPointToIndex = ...
                    objHdl.m_IndexToPhysicalPoint.^(-1);
            end
        end

        % set the spacing for the image in "mm"
        % check for 0.0 spacing values (impossible - generate error)
        % calculate the index to physical point conversion matrix
        % calculate the physical point to index conversion matrix
        function SetSpacing( objHdl, spacing )
            if spacing(:) == 0.0
                %GENERATE ERROR;
                DEBUG('ERROR - Spacing cannot be zero');
            elseif objHdl.m_Spacing ~= spacing
                objHdl.m_Spacing = spacing;
                objHdl.m_IndexToPhysicalPoint = ...
                    objHdl.m_Direction * objHdl.m_Spacing;
                objHdl.m_PhysicalPointToIndex = ...
                    objHdl.m_IndexToPhysicalPoint.^(-1);
            end
        end

        %key functions
        %---------------
        %transform indx to phy point
        % transform phys point to index
        %trans ohysical indx to physical point
        %transform physical point to continueous indx
        %cheater convenience function - get interpolated value
            %at continuous indx
            %at physical point
        %

        %---------




        % Testing '==' logic for comparison of objects
        function tf = eq(a,b) %logic for comparison
            tf = false;
            if isa(a,'MatlabITKImageObject') && isa(b,'MatlabITKImageObject')
                propsa = properties(a);
                propsb = properties(b);
                if length(propsa) == length(propsb)
                    tf = true;
                    for i=1:length(propsa)
                        %flag = COMPAREVARS ( A, B, DCP, NaNCheck, debug )
                        if ~comparevars( b.(propsb{i}), ...
                                a.(propsa{i}),1e-8,true,false )
                            tf = false;
                        end
                        if ~strcmp(propsb{i},propsa{i})
                            tf = false;
                        end
                    end
                end
            end
        end%end eq function

    end%end public member fucntions section
end

% Print Debugging material to the screen for printf style debugging
function DEBUG(input)
    if(true)
        disp(input)
    end
end

%GET OUT OF HERE PUT IN NOTES FILE

% SETTERS NO LOGIC HERE >>>>>>>>>>>>>>>>
%
%         function SetModality( objHdl, modality )
%             objHdl.m_Modality = modality;
%         end
%
%         function SetBValue( objHdl, bvalue )
%             objHdl.m_BValue = bvalue;
%         end
%
%         function SetGradientDirections( objHdl, gradientdirections )
%             objHdl.m_GradientDirections = gradientdirections;
%         end
%
%         function SetData( objHdl, data )
%             DEBUG(size(data));
%             objHdl.m_Data = data;
%             DEBUG(size(objHdl.m_Data));
%         end
% %
%         function SetLargestPossibleRegion( objHdl, region )
%             objHdl.m_LargestPossibleRegion = region;
%         end
%



%             DEBUG(   'Largest Possible Region'   )
%             DEBUG(   objHdl.GetLargestPossibleRegion()  )
%             DEBUG(   'Buffered Region'   )
%             DEBUG(   objHdl.GetBufferedRegion() )
%             DEBUG(   'Requested Region'   )
%             DEBUG(   objHdl.GetRequestedRegion() )
%             DEBUG(   'Spacing'   )
%             DEBUG(   objHdl.GetSpacing()    )
%             DEBUG(   'Origin'    )
%             DEBUG(   objHdl.GetOrigin() )
%             DEBUG(   'Direction' )
%             DEBUG(   objHdl.GetDirection()  )
%             DEBUG(   'Index to Physical Point Matix' )
%             DEBUG(   objHdl.GetIndexToPhysicalPoint()   )
%             DEBUG(   'Physical Point to Index Matix' )
%             DEBUG(   objHdl.GetPhysicalPointToIndex()   )
%             DEBUG(   'Inverse Direction' )
%             DEBUG(   objHdl.GetInverseDirection()   )
%             DEBUG(   objHdl.GetBValue()   )
%             DEBUG(   objHdl.GetModality()   )
%             DEBUG(   objHdl.GetGradientDirections()   )
%             DEBUG(   size(objHdl.m_Data)   )

%GETTERS NO LOGIC HERE>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%         function data = GetData( objHdl )
%             DEBUG( objHdl.m_Data );
%             data = objHdl.m_Data;
%             DEBUG( data );
%         end
%
%         function modality = GetModality( objHdl )
%             modality = objHdl.m_Modality;
%         end
%
%         function bvalue = GetBValue( objHdl )
%             bvalue = objHdl.m_BValue;
%         end
%
%         function gradientdirections = GetGradientDirections( objHdl )
%            gradientdirections = objHdl.m_GradientDirections;
%         end
%
%         function spacing = GetSpacing( objHdl )
%             spacing = objHdl.m_Spacing;
%         end
%
%         function origin = GetOrigin( objHdl )
%             origin = objHdl.m_Origin;
%         end
%
%         function direction = GetDirection( objHdl )
%             direction = objHdl.m_Direction;
%         end
%
%         function inverseDirection = GetInverseDirection( objHdl )
%             inverseDirection = objHdl.m_InverseDirection;
%         end
%
%         function region = GetLargestPossibleRegion( objHdl )
%             region = objHdl.m_LargestPossibleRegion;
%         end
%
%         function region = GetBufferedRegion( objHdl )
%             region = objHdl.m_BufferedRegion;
%         end
%
%         function region = GetRequestedRegion( objHdl )
%             region = objHdl.m_RequestedRegion;
%         end
%
%         function IdxToPhyPtMx = GetIndexToPhysicalPoint( objHdl )
%             IdxToPhyPtMx = objHdl.m_IndexToPhysicalPoint;
%         end
%
%         function PhyPtToIdxMx = GetPhysicalPointToIndex( objHdl )
%             PhyPtToIdxMx = objHdl.m_PhysicalPointToIndex;
%         end
%---------------------------------------------------------------
% PROBABLY UNEEDED
%         function basic = ToMatlabStruct( objHdl ) %constructor
%             basic.space = objHdl.GetSpacing();
%             basic.spaceorigin = objHdl.GetOrigin();
%             basic.spacedirections  = objHdl.GetDirection();
%             %ITKImage.SetLargestPossibleRegion(  basic.size );
%             basic.data = objHdl.GetData();
%             %basic.space
%             %basic.spacedirections
%             %basic.centerings
%             %basic.kinds
%             %basic.spaceunits
%             %basic.measurementframe
%             if(objHdl.GetBValue() > 0)
%                 basic.modality = objHdl.GetModality();
%                 basic.bvalue = objHdl.SetBValue();
%                 basic.gradientdirections = objHdl.SetGradientDirections();
%             end
%         end





%         function objHdl = FromNRRDStruct( basic , saveNRRD) %constructor
%             if nargin == 1 || saveNrrd
%                 nrrdNames = fieldnames(basic);
%                 for i = 1 : length(nrrdNames)
%                     objHdl.(nrrdNames{i}) = basic.(nrrdNames{i})
%                 end
%             end
%
%             objHdl.SetOrigin( basic.spaceorigin );
%
% %             %objHdl.3 = basic.space;
% %             objHdl.SetSpacing( basic.space );
% %             objHdl.SetOrigin( basic.spaceorigin );
% %             objHdl.SetDirection( basic.spacedirections );
% %             %ITKImage.SetLargestPossibleRegion(  basic.size );
% %             objHdl.SetData(basic.data);
% %             objHdl.m_Space = basic.space;
% %             %objHdl.m_SpaceDirections = basic.spacedirections;
% %             %objHdl.m_Centerings = basic.centerings;
% %             %objHdl.m_Kinds = basic.kinds;
% %             %objHdl.m_SpaceUnits = basic.spaceunits;
% %             %objHdl.m_SpaceDefinition = basic.spacedefinition;
% %             %objHdl.m_SpaceOrigin = basic.spaceorigin;
% % %             objHdl.m_MeasurementFrame = basic.measurementframe;
% %             %if(basic.bvalue > 0)
% %                 %objHdl.SetModality(basic.modality);
% %                 %objHdl.SetBValue(basic.bvalue);
% %                 %objHdl.SetGradientDirections(basic.gradientdirections);
% %             %end
% %             DEBUG('end of constructuor');
% %             objHdl.PrintSelf();
% %             DEBUG('end of constructor');
%         end
