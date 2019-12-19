classdef my_LinOpSelectorObj < LinOpSelector
    % Patch Selector linear operator which extracts a plane from a given
    % hypercube
    %
    % :param sz:  size of \\(\\mathrm{x}\\) = [nb_y, nb_x, nb_1, nb_2].
    % :param id_1: array containing the indexes kept along the third
    % dimension     
    % :param id_2: array containing the indexes kept along the fourth
    % dimension   
    %
    % All attributes of parent class :class:`LinOpSelector` are inherited. 
    %
    % **Example** S=my_LinOpSelectorObj(sz, id_1, id_2)
    %
    % See also :class:`LinOp`, :class:`LinOpSelector`,
    % :class:`LinOpDownsampling`.
    
    %%    Copyright (C) 2017 
    %     Modified: 07/10/2018 (mm/dd/yyyy)
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
        id_1 ;
        id_2 ;
        nb_p = 1 ;
    end
    
    %% Constructor
    methods
        function this = my_LinOpSelectorObj(sz, id_1, id_2)
            
            this.name ='LinOp my_LinOpSelectorObj';	
            assert(cmpSize(size(id_1),size(id_2)), ...
                'Parameters id_1 and id_2 must have the same size');
            this.sizein=sz;
            this.nb_p = length(id_1) ;
            this.sizeout = [this.sizein(1:2), this.nb_p];
            
            this.id_1 = id_1 ;
            this.id_2 = id_2 ;
        end
    end
    
    %% Core Methods containing implementations (Protected)
	methods (Access = protected)		
        function y = apply_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.           
            y = zeros(this.sizeout) ;
            for pl = 1:this.nb_p
                y(:, :, pl) = x(:, :, this.id_1(pl), this.id_2(pl)) ;
            end
            
        end        
        function y = applyAdjoint_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            y = zeros(this.sizein);
            for pl = 1:this.nb_p
                y(:, :, this.id_1(pl), this.id_2(pl)) = x(:, :, pl) ;
            end
        end
        function y = applyHtH_(this,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            y = zeros(this.sizein);
            for pl = 1:this.nb_p
                y(:, :, this.id_1(pl), this.id_2(pl)) = ...
                    x(:, :, this.id_1(pl), this.id_2(pl)) ;
            end       
        end
        function y = applyHHt_(~,x)
            % Reimplemented from parent class :class:`LinOpSelector`.  
            y = x;
        end
        function M = makeHtH_(this)
            % Reimplemented from parent class :class:`LinOpSelector`. 
            w=zeros(this.sizein);
            for pl = 1:this.nb_p
                w(:, :, this.id_1(pl), this.id_2(pl)) = 1 ;
            end
            M=LinOpDiag([],w);
        end
    end
end