function [Hcl,sort_all_gm,sort_all_pm,mu_results,all_lp] ...
    = do_RSanal(A_cl,B_cl,C_cl,D_cl,u_names,y_names,ww,do_prt,do_plt,w_scl,freq_units)
% outputs of interest: 
% Hcl,all_lp
% w_mu_u,w_mu_y,w_mu_uy,
% max_ubnds,max_ybnds,max_uybnds,
% w_mu_max,sort_all_gm,sort_all_pm,
% igm_gmin1,wgm_gmin1,gm_gmin1,igm_gmin2,wgm_gmin2,gm_gmin2,
% ipm_pmin1,wpm_pmin1,pm_pmin1,ipm_pmin2,wpm_pmin2,pm_pmin2

if( 0 )
    % an example but
    % u_names and y_names should be set outside
    % ww should be set outside
    % do_prt should be set outside
    % do_plt should be set outside
    % w_scl should be set outside
    % freq_units should be set outside
    u_names = {'da','de','dr','throt' } ; %#ok<UNRCH>
    y_names = { 'p','q','r','phi','theta','psi','Nd','Ed','hd','N','E','h','axa','aya','aza','V','alpha','beta' } ;
    ww = logspace(-1,1,200) ;
    do_prt = 1 ;
    do_plt = 1 ;
    w_scl = 1 ;
    freq_units = 'rad/sec' ;
end

% always subplot(2,1,:) gain for u and y 
plot_phase = 0 ; % subplot(2,1,:) phase for u and y loops
plot_gain_and_phase = 0 ; % subplot(2,1,:) gain and phase in same figure

% select upper and lower limits for gain plots (to focus on gain crossover)
plt_mag_min = 0.01 ;
plt_mag_max = 100 ;

% threshold to ignore loops with small gain
abslp_threshold = 1.e-3 ;

% to reduce number of figures when there are many actuators and sensors
% all u in one figure, y(1:6), y(7:12), y(13:18) in separate figures
use_several_plot_pages = 0 ;
    
num_u = length(u_names) ;
num_y = length(y_names) ;
num_uy = num_u + num_y ;

uy_names = horzcat(u_names,y_names) ;

clsys = ss(A_cl,B_cl,C_cl,D_cl) ;

% do the Frequency response analysis
Hcl = freqresp(clsys,ww) ;

num_ww = length(ww) ;
Hu = NaN(num_u,num_ww) ;
Hy = NaN(num_y,num_ww) ;
Huy = NaN(num_uy,num_ww) ;

i_u = 1:num_u ;
i_y = num_u+1:num_uy ;
i_uy = [ i_u i_y ] ;
for i=1:num_ww
    Hu(:,i) = svd(squeeze(Hcl(i_u,i_u,i))) ;
    Hy(:,i) = svd(squeeze(Hcl(i_y,i_y,i))) ;
    Huy(:,i) = svd(squeeze(Hcl(i_uy,i_uy,i))) ;
end

disp(' ')
[Hu_bnds,MUINFO] = mussv(Hcl(i_u,i_u,:),ones(num_u,2)) ;   %#ok<ASGLU>
[Hy_bnds,MUINFO] = mussv(Hcl(i_y,i_y,:),ones(num_y,2)) ;   %#ok<ASGLU>
[Huy_bnds,MUINFO] = mussv(Hcl(i_uy,i_uy,:),ones(num_uy,2)) ;   %#ok<ASGLU>
% 1 x 2 x num_ww

if( num_ww > 1 )
    Muu = squeeze(Hu_bnds)' ; % num_ww x 2
    Myy = squeeze(Hy_bnds)' ;
    Muy = squeeze(Huy_bnds)' ;
else
    Muu = squeeze(Hu_bnds) ; % num_ww x 2
    Myy = squeeze(Hy_bnds) ;
    Muy = squeeze(Huy_bnds) ;
end

[max_ubnds,Iu_bnds] = max(Muu(:,1),[],1) ; %#ok<ASGLU>
[max_ybnds,Iy_bnds] = max(Myy(:,1),[],1) ; %#ok<ASGLU>
[max_uybnds,Iuy_bnds] = max(Muy(:,1),[],1) ; %#ok<ASGLU>
w_mu_u = ww(Iu_bnds) ;  %#ok<NASGU>
w_mu_y = ww(Iy_bnds) ;  %#ok<NASGU>
w_mu_uy = ww(Iuy_bnds) ;  %#ok<NASGU>

[maxMuu,indx_Muu] = max(Muu,[],1) ; % 1 x 2
[maxMyy,indx_Myy] = max(Myy,[],1) ;
[maxMuy,indx_Muy] = max(Muy,[],1) ;
w_mu_max = [ ww(indx_Muu(1))*w_scl maxMuu(1) ;  ...
    ww(indx_Myy(1))*w_scl maxMyy(1) ; ...
    ww(indx_Muy(1))*w_scl maxMuy(1) ] ;  %#ok<NASGU>
max_over_w_bnds = [ maxMuu ; maxMyy ; maxMuy ] ; % 3x2

mu_results = [ ww(indx_Muu(1))*w_scl,maxMuu(1) ww(indx_Myy(1))*w_scl,maxMyy(1) ww(indx_Muy(1))*w_scl,maxMuy(1) ] ;

if( do_prt )
    disp('--------------------------------------------------------');
    max_over_w_lower_and_upper_bnds = max_over_w_bnds(:,[2 1]) %#ok<NASGU,NOPRT>
    disp('--------------------------------------------------------');
    disp('Loop       w      max(mu)')
    disp(sprintf('%6s %8.4f %8.2f','u' ,ww(indx_Muu(1))*w_scl,maxMuu(1))) %#ok<DSPS>
    disp(sprintf('%6s %8.4f %8.2f','y' ,ww(indx_Myy(1))*w_scl,maxMyy(1))) %#ok<DSPS>
    disp(sprintf('%6s %8.4f %8.2f','uy' ,ww(indx_Muy(1))*w_scl,maxMuy(1))) %#ok<DSPS>
end

if( do_plt && num_ww > 1 )
    figure
    subplot(3,1,1) , loglog(ww*w_scl,squeeze(Hu_bnds(1,1:2,:)), ...
        [ww(1) ww(end)]*w_scl,[1 1],'k-')
    axis([ ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
    title('Structured Singular Value upper and lower bounds for multiplicative perturbations')
    ylabel([ num2str(num_u) ' Actuator'])
    subplot(3,1,2) , loglog(ww*w_scl,squeeze(Hy_bnds(1,1:2,:)), ...
        [ww(1) ww(end)]*w_scl,[1 1],'k-')
    axis([ ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
    ylabel([ num2str(num_y) ' Sensors'])
    subplot(3,1,3) , loglog(ww*w_scl,squeeze(Huy_bnds(1,1:2,:)), ...
        [ww(1) ww(end)]*w_scl,[1 1],'k-')
    axis([ ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
    ylabel([ num2str(num_uy) ' Actuator + Sensors'])
    xlabel(['Frequency (' freq_units ')'])

    figure
    loglog(ww*w_scl,squeeze(Hu_bnds(1,1,:)),'b-', ...
        ww*w_scl,squeeze(Hy_bnds(1,1,:)),'g--', ...
        ww*w_scl,squeeze(Huy_bnds(1,1,:)),'r-.', ...
        [ww(1) ww(end)]*w_scl,[1 1],'k-')
    title('Structured Singular Value upper bounds for multiplicative perturbations')
    axis([ ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
    ylabel('Upper bound for ssv')
    xlabel(['Frequency (' freq_units ')'])
    legend([num2str(num_u) ' Actuator'],[num2str(num_y) ' Sensors'],[num2str(num_uy) ' Actuator + Sensors'])
    % legend('1 actuator','1 sensor') % did get help from spratt on 2/6/2014
end

% calculate L = lp = loop transfer function for negative feedback
% starting with closed loop transfer function = M = Hcl(k,k)
% assume input and output k of closed loop correspond to mult. perturbation
% y = - L * ( u + y )
% (1+L)*y = - L * u
% y = - L/(1+L) * u
% let M = -L/(1+L)
% (1+L) * M = -L
% (1 + M ) * L + M = 0
% L = - M / ( 1 + M )

all_gm = [] ;
all_pm = [] ;
all_lp = NaN(num_ww,num_uy) ;

for idx = 1:num_uy
    mcl = squeeze(Hcl(idx,idx,:)) ;
    lp = - ( 1 + mcl ) .\ mcl ;
    if( max(abs(lp)) > abslp_threshold )
        [ gm_results,pm_results ] = gpm([lp ww'],idx) ;
        all_gm = [ all_gm ; gm_results ] ; %#ok<AGROW>
        all_pm = [ all_pm ; pm_results ] ; %#ok<AGROW>
    end
    all_lp(:,idx) = lp ;
end

if( ~isempty(all_gm) )
    all_gm = all_gm * diag([1 w_scl 1 ]) ;
end

if( ~isempty(all_pm) )
    all_pm = all_pm * diag([1 w_scl 1 ]) ;
end

name_length = -max(cellfun('length',uy_names)) ;
if( do_prt )
    disp('--------------------------------------------------------');
    disp(' ')
    disp('Loop      w        GM')
    disp(' ')
end
if( size(all_gm,1) > 0 )
    disp('--------------------------------------------------------');
    [Y,I_gm] = sort(abs(all_gm(:,3)),1,'ascend') ; %#ok<ASGLU>
    sort_all_gm = all_gm(I_gm,:) ;
    if( do_prt )    
        for kkk = I_gm'
            disp(sprintf('%*s %8.4f %8.2f',name_length,uy_names{all_gm(kkk,1)},all_gm(kkk,2:3))) %#ok<DSPS>
        end        
    end
    % store index, Frequency, smallest 2 margins
    igm_gmin1 = all_gm(I_gm(1),1) ; %#ok<NASGU>
    wgm_gmin1 = all_gm(I_gm(1),2) ; %#ok<NASGU>
    gm_gmin1 = all_gm(I_gm(1),3) ; %#ok<NASGU>
    if( size(all_gm,1) > 1 )
        igm_gmin2 = all_gm(I_gm(2),1) ; %#ok<NASGU>
        wgm_gmin2 = all_gm(I_gm(2),2) ; %#ok<NASGU>
        gm_gmin2 = all_gm(I_gm(2),3) ; %#ok<NASGU>
    else
        igm_gmin2 = NaN ; %#ok<NASGU>
        wgm_gmin2 = NaN ; %#ok<NASGU>
        gm_gmin2 = Inf ; %#ok<NASGU>
    end
else
    disp('--------------------------------------------------------');
    disp('No gm found')
    sort_all_gm = [ ] ;
end

if( do_prt )
    disp('--------------------------------------------------------');
    disp(' ')
    disp('Loop      w        PM')
    disp(' ')
end
if( size(all_pm,1) > 0 )
    disp('--------------------------------------------------------');
    [Y,I_pm] = sort(abs(all_pm(:,3)),1,'ascend') ; %#ok<ASGLU>
    sort_all_pm = all_pm(I_pm,:) ;
    if( do_prt )
        for kkk = I_pm'
            disp(sprintf('%*s %8.4f %8.2f',name_length,uy_names{all_pm(kkk,1)},all_pm(kkk,2:3))) %#ok<DSPS>
        end
    end
    ipm_pmin1 = all_pm(I_pm(1),1) ; %#ok<NASGU>
    wpm_pmin1 = all_pm(I_pm(1),2) ; %#ok<NASGU>
    pm_pmin1 = all_pm(I_pm(1),3) ; %#ok<NASGU>
    if( size(all_pm,1) > 1 )
        ipm_pmin2 = all_pm(I_pm(2),1) ; %#ok<NASGU>
        wpm_pmin2 = all_pm(I_pm(2),2) ; %#ok<NASGU>
        pm_pmin2 = all_pm(I_pm(2),3) ; %#ok<NASGU>
    else
        ipm_pmin2 = NaN ; %#ok<NASGU>
        wpm_pmin2 = NaN ; %#ok<NASGU>
        pm_pmin2 = NaN ; %#ok<NASGU>
    end
else
    disp('--------------------------------------------------------');
    disp('no pm found')
    sort_all_pm = [ ] ;
end

if( do_plt && num_ww > 1 )
    if( use_several_plot_pages == 1 )
        disp('--------------------------------------------------------');
        figure
        loglog(ww*w_scl,abs(all_lp(:,i_u)),[ww(1) ww(end)]*w_scl,[1 1],'k-')
        title('Broken loop single-loop-at-a-time analysis')
        axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
        ylabel('Loop gain')
        legend(u_names)

        figure
        yidx = 1:6 ;
        pltidx = num_u + yidx ;
        loglog(ww*w_scl,abs(all_lp(:,pltidx)),[ww(1) ww(end)]*w_scl,[1 1],'k-')
        title('Broken loop single-loop-at-a-time analysis')
        axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
        ylabel('Loop gain')
        xlabel(['Frequency (' freq_units ')'])
        legend(y_names(yidx))
        grid on
        
        yidx = 7:12 ;
        if( num_y >= yidx(end) )
            figure
            pltidx = num_u + yidx ;
            loglog(ww*w_scl,abs(all_lp(:,pltidx)),[ww(1) ww(end)]*w_scl,[1 1],'k-')
            title('Broken loop single-loop-at-a-time analysis')
            axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
            ylabel('Loop gain')
            xlabel(['Frequency (' freq_units ')'])
            legend(y_names(yidx))
            grid on
        end
        
        yidx = 13:18 ;
        if( num_y >= yidx(end) )
            figure
            pltidx = num_u + yidx ;
            loglog(ww*w_scl,abs(all_lp(:,pltidx)),[ww(1) ww(end)]*w_scl,[1 1],'k-')
            title('Broken loop single-loop-at-a-time analysis')
            axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
            ylabel('Loop gain')
            xlabel(['Frequency (' freq_units ')'])
            legend(y_names(yidx))
            grid on
        end
    else
        figure
        % eliminate plots of loop gains with gain less than what is plotted
        plt_uidx = [ ] ;
        for kk = 1:length(i_u)
            if( max(abs(all_lp(:,i_u(kk)))) > plt_mag_min )
                plt_uidx = [ plt_uidx kk ] ; %#ok<AGROW>
            end
        end
        plt_yidx = [ ] ;
        for kk = 1:length(i_y)
            if( max(abs(all_lp(:,i_y(kk)))) > plt_mag_min )
                plt_yidx = [ plt_yidx kk ] ; %#ok<AGROW>
            end
        end
        subplot(2,1,1)
        if( ~isempty(plt_uidx) )
            loglog(ww*w_scl,abs(all_lp(:,i_u(plt_uidx))), ...
                [ww(1) ww(end)]*w_scl,[1 1],'k-')
            title('Broken loop single-loop-at-a-time analysis')
            axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
            ylabel('Loop gain')
            legend(u_names(plt_uidx))
        end
        subplot(2,1,2)
        if( ~isempty(plt_uidx) )
            loglog(ww*w_scl,abs(all_lp(:,i_y(plt_yidx))), ...
                [ww(1) ww(end)]*w_scl,[1 1],'k-')
            axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
            ylabel('Loop gain')
            xlabel(['Frequency (' freq_units ')'])
            legend(y_names(plt_yidx))
        end
        grid on
        
        if( plot_phase == 1 )
            figure
            subplot(2,1,1) , semilogx(ww*w_scl,angle(all_lp(:,i_u))*180/pi, ...
                [ww(1) ww(end)]*w_scl,[-180 -180],'k-')
            title('Broken loop single-loop-at-a-time analysis')
            %             axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
            ylabel('Loop phase')
            legend(u_names)
            subplot(2,1,2) , semilogx(ww*w_scl,angle(all_lp(:,i_y))*180/pi, ...
                [ww(1) ww(end)]*w_scl,[-180 -180],'k-')
            %             axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
            ylabel('Loop phase (deg)')
            xlabel(['Frequency (' freq_units ')'])
            legend(y_names)
            grid on
        end
        
        if( plot_gain_and_phase == 1 )
            for i_lp = [ i_u i_y ]
                if( max(abs(all_lp(:,i_lp))) >= 1.e-5)
                    figure
                    subplot(2,1,1) , loglog(ww*w_scl,abs(all_lp(:,i_lp)), ...
                        [ww(1) ww(end)]*w_scl,[1 1],'k-')
                    idx_u_intersect = intersect( i_lp, i_u ) ;
                    idx_y_intersect = intersect( i_lp, i_y ) ;
                    if( ~isempty(idx_u_intersect) )
                        idx_u = find( i_lp == i_u ) ;
                        %                     title(['broken loop single-loop-at-a-time analysis for u loop ' num2str(idx_u)])
                        title(['Broken loop single-loop-at-a-time analysis for ' ...
                            u_names{idx_u} ' loop']) %#ok<FNDSB>
                    elseif( ~isempty(idx_y_intersect) )
                        idx_y = find( i_lp == i_y ) ;
                        %                     title(['broken loop single-loop-at-a-time analysis for y loop ' num2str(idx_y)])
                        title(['Broken loop single-loop-at-a-time analysis for ' ...
                            y_names{idx_y} ' loop']) %#ok<FNDSB>
                    end
                    axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
                    ylabel('Loop gain')
                    subplot(2,1,2) , semilogx(ww*w_scl,angle(all_lp(:,i_lp))*180/pi, ...
                        [ww(1) ww(end)]*w_scl,[-180 -180],'k-')
                    %             axis([ww(1)*w_scl ww(end)*w_scl plt_mag_min plt_mag_max ])
                    ylabel('Loop phase (deg)')
                    xlabel(['Frequency (' freq_units ')'])
                    grid on
                end
            end
        end
    end
end
