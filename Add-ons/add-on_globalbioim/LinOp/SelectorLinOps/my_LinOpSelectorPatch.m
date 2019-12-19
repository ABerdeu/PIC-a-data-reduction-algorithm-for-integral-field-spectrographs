classdef my_LinOpSelectorPatch < LinOpSelector
    % Patch Selector linear operator which extracts a patch from the given
    % vector
    % $$\\mathrm{H} : \\mathrm{x} \\mapsto \\mathrm{x}_{[i_{min}:i_{max}]}$$
    % where \\( i_{min} \\) and \\( i_{max} \\) are indexes corresponding
    % to the first and last elements of the patch.
    %
    % :param sz:  size of \\(\\mathrm{x}\\) on which the :class:`LinOpDownsample` applies.
    % :param idxmin: array containing the first kept index in each direction         
    % :param idxmax: array containing the last kept index in each direction
    % :param flag_squeeze: boolean to squeeze (true) or not (false)
    % singleton dimensions (default false)
    %
    % All attributes of parent class :class:`LinOpSelector` are inherited. 
    %
    % **Example** S=my_LinOpSelectorPatch(sz, idxmin, idxmax)
    %
    % **Example** S=my_LinOpSelectorPatch(sz, idxmin, idxmax, flag_squeeze)
    %
    % See also :class:`LinOp`, :class:`LinOpSelector`,
    % :class:`LinOpDownsampling`.
    
    %%    Copyright (C) 2017 
    %     E. Soubies  emmanuel.soubies@epfl.ch
    %     Modified: 06/05/2018 (mm/dd/yyyy)
    %     Anthony Berdeu (Laboratoire Hubert Curien)
    %
    %     This program is free software: you can redistribute it and/or modify
    %     it under the terms of the GNU General Public License as published by
    %     the Free Software Foundation, either version 3 of the License, or
    %     (at your option) any later version.
    %
    %     This program is distributed in the hope that it will be useful,
    %     but WITHOUT ANY WARRANTY; without even the implied warranty of
    %     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    %     GNU General Public License for more details.
    %
    %     You should have received a copy of the GNU General Public License
    %     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    properties (SetAccess = protected,GetAccess = public)
        flag_squeeze ; % boolean to squeeze (true) or not (false)
            % singleton dimensions
    end
    
    %% Constructor
    methods
        function this = my_LinOpSelectorPatch(sz, idxmin, idxmax, ...
                flag_squeeze)
            
            % flag_squeeze
            if nargin<4
                flag_squeeze = false ;
            elseif isempty(flag_squeeze)
                flag_squeeze = false ;
            end
            this.flag_squeeze = flag_squeeze ;
            
            this.name ='LinOp my_SelectorPatch';	
            assert(cmpSize(size(sz),size(idxmin)) && ...
                cmpSize(size(sz),size(idxmax)), ...
                'Parameters sz idxmin and idxmax must have the same size');
            assert(~any(idxmin<=0) && ~any(idxmin>sz), ...
                'idxmin out of sz range');
            assert(~any(idxmax<=0) && ~any(idxmax>sz), ...
                'idxmax out of sz range');
            assert(~any(idxmin>idxmax), ...
                'idxmin must be smaller than idxmax in each dimension');
            this.sizeout=idxmax-idxmin+1;
            if this.sizeout(end)==1, this.sizeout=this.sizeout(1:end-1);end;
            if flag_squeeze
                this.sizeout = this.sizeout(~(this.sizeout==1)) ;
            end
            while length(this.sizeout)<2 % Taking into account scalars and
                % vectors
                this.sizeout = [this.sizeout, 1] ;
            end
            this.sizein=sz;
            this.sel=cell(length(sz),1);
            for ii=1:length(sz)
                this.sel{ii}=idxmin(ii):idxmax(ii);
            end
        end
    end
    
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)		
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.           
            y =x(this.sel{:});
            if this.flag_squeeze
                y = squeeze(y) ;
                if size(y, 1)==1 % y is row-vector. Shaping it into a
                        % column-vector 
                    y = y.' ;
                end
            end
        end        
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            y = zeros(this.sizein);
            y(this.sel{:}) = x;
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            y = zeros(this.sizein);
            y(this.sel{:}) = x(this.sel{:});            
        end
        function y = applyHHt_(~,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            y = x;
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOpSelector`. 
            w=zeros(this.sizein);w(this.sel{:})=1;
            M=LinOpDiag([],w);
        end
    end
end