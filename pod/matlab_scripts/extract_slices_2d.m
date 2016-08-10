
[data,lr1,elmap,time,istep,fields,emode,wdsz,etag,header,status] = readnek('./pipe_mode_sym0.f00007');
if status < 0
    error('Could not read file');
end

%z_val=[-2 -1 -.5:0.25:4 5 8 12]; %! corresponds to upstream section, exptected to be <= 0
z_val=[-2:22];

nFpp = 176;
bent_phi=90*pi/180; % bent section
bent_radius=5./3.; %1.66666666666666666666666666667, ! bent radius (defined from pipe centerline)
circumf = bent_radius*bent_phi;
c_ = cos(bent_phi);
s_ = sin(bent_phi);

data2d = zeros(nFpp,lr1(1)*lr1(2),6);

for m=1:(size(data,1)/nFpp)
    for gll = 1:lr1(1)-1
        
        for ex=1:nFpp
            e  = (m-1)*nFpp+ex;
            
            %   - data:   nek5000 data ordered as (iel,inode,[x|y|(z)|u|v|(w)|p|T|s_i])
            
            z_unbent = data(e,(gll-1)*lr1(1)*lr1(2)+1,3);
            c = 1;
            s = 0;
            
            if (abs(bent_phi) > 1e-6)
                if (z_unbent > 0)
                    angle = atan2(z_unbent,data(e,(gll-1)*lr1(1)*lr1(2)+1,1));
                    if (angle <= bent_phi)
                        z_unbent = bent_radius*angle;
                        c = cos(-angle);
                        s = sin(-angle);
                    else
                        z_unbent = c_*z_unbent/s_^2 - data(e,(gll-1)*lr1(1)*lr1(2)+1,1)/s_;
                        z_unbent = z_unbent/(1+c_^2/s_^2) + circumf;
                        c = cos(-bent_phi);
                        s = sin(-bent_phi);
                    end
                end
            end
            
            % numerical precission
            for k=1:(lr1(1)*lr1(2))
                
                data2d(ex,k,1) = c*data(e,(gll-1)*lr1(1)*lr1(2)+k,1)-s*data(e,(gll-1)*lr1(1)*lr1(2)+k,3); %x
                data2d(ex,k,2) = data(e,(gll-1)*lr1(1)*lr1(2)+k,2); %y
                
                
                data2d(ex,k,3) = c*data(e,(gll-1)*lr1(1)*lr1(2)+k,4)-s*data(e,(gll-1)*lr1(1)*lr1(2)+k,6); %u
                
                data2d(ex,k,4) = data(e,(gll-1)*lr1(1)*lr1(2)+k,5); %v
                data2d(ex,k,5) = s*data(e,(gll-1)*lr1(1)*lr1(2)+k,4)-c*data(e,(gll-1)*lr1(1)*lr1(2)+k,6);
                
                data2d(ex,k,6) = 0; %T
            end
            
        end
        data2d(:,:,1) = data2d(:,:,1) - .5*(max(max(data2d(:,:,1)))+min(min(data2d(:,:,1))));

        lr1(3)=1;
        elmap = int32([1:nFpp]);
        writenek(['mode6/slab_mode_2d0.f' num2str((m-1)*(lr1(1)-1)+gll,'%05d')],data2d,lr1,elmap,z_unbent,m,fields,emode,wdsz,etag);
    end
end

