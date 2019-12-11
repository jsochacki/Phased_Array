function [FileName, sp, freq_GHz] = ReadTouchstone()
% function [sp,freq]=ReadTouchsone()
% Reads the touchstone .snp file
% Returns S parameters in complex 3D matrix in sp(n,m,freq)
% Returns frequency in freq in [GHz]
% Value of n in .snp must be the number of ports
% Normalization factor is not returned in this script (the number after
% R in the header), adding z0 in the return value should give this back

% Touchstone format is based on the version 1.1
% http://www.eda.org/pub/ibis/connector/touchstone_spec11.pdf

% base on the script found on the edaboard.com and also N.K.
% http://www.edaboard.com/ftopic328289.html
% Ryosuke Ito on 12/02/2009



verbose=1 ; % flag for verbose mode, if other than 0, shows the progress

[FileName,PathName]=uigetfile('*.s*p','Select valid touchstone file');
% if cancel is pressed, returns zeros
if FileName==0
    sp=0 ; freq_GHz=0 ;
    return
end

[t,r]=strtok(FileName,'.'); % separating before and after the first '.'
[t,r]=strtok(r,'.'); % extracting the first token after the first '.'
s1=sscanf(t,'%c%d%c'); % extracting the integer in between the characters from extension
port_num=s1(2);  % getting the n in 'snp'

fname=[PathName, FileName];
fid=fopen(fname,'r');

% line_entries is the number of complex data per line
if(port_num>4)
    line_entries=4;
    positions_per_last_line=rem(port_num,4);
elseif(port_num==2)
    line_entries=4;
else
    line_entries=port_num;
end

%%%%%%%%%%%%%%%%%%%
mode=0; % data format mode 'RI' or 'MA'
unit=1; % frequency unit, 1 for GHz, 1e-3 for MHz
j=sqrt(-1); % defining j


%% process the preamble.
buffer='NaN' ;

while( isnan(str2double(strtok(buffer))) ) % checking if the first token is not a number
    
    % read one line from the file, change to upper case, and save to buffer
    buffer=upper(fgetl(fid));  

    % detecting blank line and discard
    if(strcmp(strtok(buffer),''))
        buffer=fgetl(fid);
        continue;
    end;
    
    %    if the first non space character is !    
    if(strncmpi(strjust(buffer,'left'),'!',1) )
        continue;
    end;
    
    % if the first non space character is #
    if(strncmpi(strjust(buffer,'left'),'#',1) )
        [token, buffer]=strtok(buffer); % dropping the # and forward the buffer
        [token, buffer]=strtok(buffer); % reading the next token
        %
        if (strcmp(token, 'HZ'))
            unit=1.0E-9;
        elseif (strcmp(token, 'MHZ'))
            unit=1.0E-3;
        elseif (strcmp(token, 'GHZ'))
            unit=1.0;
        end;
        %
        [token, buffer]=strtok(buffer); % dropping 'S' 
        [token, buffer]=strtok(buffer);
        %
        if (strcmp(token, 'RI'))
            mode=1;
         elseif (strcmp(token, 'MA'))
            mode=0;

        end;
        [token, buffer]=strtok(buffer); % dropping 'R'
        [token, buffer]=strtok(buffer);
        z0=str2num(token) ; % normalization value saved into z0
    end
    
end

% for verbose mode
switch verbose
    case 0 
    otherwise
        switch mode
            case 1
                disp('Reading touchstone file in RI mode.') ;
            case 0
                disp('Reading touchstone file in MA mode.') ;
        end
        
        switch unit
            case 1
                disp('Frequency unit is GHz.') ;
            case 1e-3
                disp('Frequency unit is MHz.') ;
            case 1e-9
                disp('Frequency unit is Hz.') ;
        end
        
        fprintf('Z0 = %f\n\n',z0) ;
end


%% Reading the data.
n=1;
while(strcmp(strtok(buffer),''))
   buffer=fgetl(fid);
end
while(~feof(fid))    
    [data,buffer]=strtok(buffer) ;
    [data_x, count]=sscanf(data, '%f', 1);
    freq_GHz(n,1)=data_x*unit;
    
    l=1;
    while (l<port_num+1)
        m=1;
        while (m<port_num+1)
            [data, buffer]=strtok(buffer);
            [data_x, count]=sscanf(data, '%f', 1);
            [data, buffer]=strtok(buffer);
            [data_y, count]=sscanf(data, '%f', 1);
            %
            if(count==0)
                buffer=fgetl(fid);
                continue;
            end;
            %
            if(mode==1)
                % when number of ports is 2, changes data order
                % see Touchstone(TM) specification
                if(port_num==2) 
                    sp(m,l,n)=data_x+j*data_y;
                else
                    sp(l,m,n)=data_x+j*data_y;
                end;
            else
                if(port_num==2)
                    sp(m,l,n)=data_x*exp(j*data_y*pi/180);
                else
                    sp(l,m,n)=data_x*exp(j*data_y*pi/180);
                end;
            end;
            m=m+1;
        end;
        l=l+1;
    end;
    n=n+1;
    
    % detecting blank lines b/w frequencies (often produced by HFSS)
    while(strcmp(strtok(buffer),''))
        buffer=fgetl(fid);
    end
end

fclose(fid);

fprintf('File read from %s\n',FileName) ;