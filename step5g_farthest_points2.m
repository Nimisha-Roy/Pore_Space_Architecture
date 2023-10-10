function [im_new,v_new]=step5g_farthest_points2(dsp,n,v,pp,c,im,debug)
%This function finds farthest points in second half of closed boundary for dominant points in step 5g
% ------
% Author: Nimisha Roy
% e-mail: nroy9@gatech.edu
% ------
%------------------------------------------------------------------------------START CODE------------------------------------------------------------------------------
dsp1=c;sd=c;p=1;p1=1;im_new=[];v_new=[];xh=[];yh=[];indexh=[];dist_p=[];dsps_store=[];
while dsp>n
        for i= 1: (size(im,1)-1)
            for j= im(i) : im(i+1)
                pt = [pp(j,2),pp(j,1),0];
                a = v(i,:)-v(i+1,:);%Vector
                b = pt-v(i+1,:);%Vector
                dist_p = [dist_p; pp(j,2) pp(j,1) (norm(cross(a,b)) / norm(a)) j];
            end
            [dsps,indh] = max(dist_p(:,3));
            dsps_store=[dsps_store;dsps];
            xh= [xh; dist_p(indh,1)];
            yh= [yh; dist_p(indh,2)];
            indexh= [indexh; dist_p(indh,4)];
            if debug==1
                disp(dist_p(indh,1));disp(dist_p(indh,2));
                h= impoint(gca, dist_p(indh,1),dist_p(indh,2));setColor(h,'m');
            end
            dist_p=[];
            if i== (size(im,1)-1)
                dsp=max(dsps_store);
                if debug ==1
                disp(['dsp is',num2str(dsp)]);
                end
                dsps_store=[];
            end
        end
        
        for l= 1: ((size(im,1))+(size(im,1)-1))
            if mod(l,2) == 1
                im_new(l,:)= im(p,:);v_new(l,:) = v(p,:);
                p=p+1;% to check how many times the loop is executed
            else
                im_new(l,:)= indexh(p1,:);
                v_new(l,1) = xh(p1,1);v_new(l,2) = yh(p1,1);
                v_new(l,3) = 0;
                p1=p1+1;% to check how many times the loop is executed
            end
        end
        im= im_new;v= v_new;xh=[];yh=[];indexh=[];p=1;p1=1;
        if dsp> n
            im_new=[];v_new=[];
        else
            for i= 1: (size(im,1)-1)
                for j= im(i) : im(i+1)
                    pt = [pp(j,2),pp(j,1),0];
                    a = v(i,:)-v(i+1,:);%Vector
                    b = pt-v(i+1,:);%Vector
                    dist_p = [dist_p; pp(j,2) pp(j,1) (norm(cross(a,b)) / norm(a)) j];
                end
                [dsp1,indh] = max(dist_p(:,3));
                if dsp1 > n
                    xh= [xh; dist_p(indh,1)];yh= [yh; dist_p(indh,2)];
                    indexh= [indexh; dist_p(indh,4)];
                else
                    xh= [xh; 0];yh= [yh;0];indexh= [indexh; 0];
                end
                dist_p=[];
            end
            for l= 1: ((size(im,1))+(size(im,1)-1))
                if mod(l,2) == 1
                    im_new(l,:)= im(p,:);v_new(l,:) = v(p,:);
                    p=p+1;% to check how many times the loop is executed
                else
                    im_new(l,:)= indexh(p1,:);
                    v_new(l,1) = xh(p1,1);v_new(l,2) = yh(p1,1);v_new(l,3) = 0;
                    p1=p1+1;% to check how many times the loop is executed
                end
            end
        end
    end
    v_new1=[];im_new1=[];
    for i= 1: size(v_new,1)
        if im_new(i,1) ~=0
            v_new1=[v_new1; v_new(i,:)];im_new1=[im_new1; im_new(i,:)];
        end
    end
    v_new= v_new1;im_new= im_new1;