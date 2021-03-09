classdef DTG<handle
    
    properties
        theta_0;
        d_theta_0;
        theta_end;
        lynxmotion;
        Object;
        Sigma;
        T_begin;
        T_end;
        P_tao;
        Point_number;
        global_point_number;
        global_sigma;
        global_t_begin;
        Threshold;
        Alpha;
        N;
        parameter;
        global_path;
        global_velo;
        global_carpoints;
        global_distance;
        time;
        state;
        init_path;
        init_velo;
        init_carpoints;
        init_para;
        init_distance;
        

    end
    methods
        function dtg = DTG(theta_0, d_theta_0, theta_end, varargin)
        
            dtg.theta_0 = theta_0;
            dtg.d_theta_0 = d_theta_0;
            dtg.theta_end = theta_end;

            definedOrDefault = @(name,default) ...
                         definedOrDefault_long(name,default,varargin);

            a1 = definedOrDefault('a1',0);
            alpha1 = definedOrDefault('alhpa1',pi/2);
            d1 = definedOrDefault('d1',0.086);

            a2 = definedOrDefault('a2',0.146);
            alpha2 = definedOrDefault('alhpa2',0);
            d2 = definedOrDefault('d2',0);

            a3 = definedOrDefault('a3',0.187);
            alpha3 = definedOrDefault('alhpa3',0);
            d3 = definedOrDefault('d3',0);
            
            dtg.N = definedOrDefault('N', 10);
            dtg.Threshold = definedOrDefault('Threshold', 0.08);
            dtg.P_tao = definedOrDefault('P_tao', 0);
            
            dtg.T_begin = definedOrDefault('T_begin',0);
            dtg.T_end = definedOrDefault('T_end',2);
            dtg.Alpha = definedOrDefault('Alpha', 200);
            dtg.Sigma = dtg.T_end - dtg.T_begin;
            
            dtg.Point_number = definedOrDefault('Point_number', 200);
            dtg.global_t_begin = definedOrDefault('global_t_begin',0);
            dtg.global_point_number = definedOrDefault('global_point_number', 200);
            dtg.global_sigma = dtg.T_end - dtg.global_t_begin;
            dtg.parameter = zeros(1,9);
            dtg.Object = definedOrDefault('Object', zeros(1,3));
            dtg.time = 0;
            dtg.global_distance = zeros(1,2);
            dtg.state = 'initialization';


            L(1) = Link('d',d1,'a', a1,'alpha', alpha1 );
            L(2) = Link('d',  d2, 'a',   a2  ,'alpha',  alpha2);
            L(3) = Link('d', d3  ,'a',  a3  ,'alpha', alpha3 );

            L(1).m = 0.3;
            L(2).m = 0.2;
            L(3).m = 0.5;
            L(1).qlim = [-pi/2, pi/2];
            L(2).qlim = [0, pi];
            L(3).qlim = [-pi, 0];
            dtg.lynxmotion = SerialLink(L, 'name', ' Lynxmotion ');
            
        end
        
        
        
        function [E_at] = objectFunction(dtg, a)

            a1 = a(1);
            a2 = a(2);
            a3 = a(3);
            a6 = a(4);
            a7 = a(5);
            a8 = a(6) ;
            a11 = a(7);
            a12 = a(8);
            a13 = a(9);

            a5 = dtg.theta_0(1);
            a10 = dtg.theta_0(2);
            a15 = dtg.theta_0(3);


            a4 = dtg.d_theta_0(1);
            a9 = dtg.d_theta_0(2);
            a14 = dtg.d_theta_0(3);

            X_obs = dtg.Object(dtg.time+1,1);
            Y_obs = dtg.Object(dtg.time+1,2);
            Z_obs = dtg.Object(dtg.time+1,3);
            
            t_begin = dtg.T_begin;
            p_tao = dtg.P_tao;
            sigma = dtg.Sigma;
   
            
            d1 = ((cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))^2*(43*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)) - 500*Z_obs*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)) + 500*X_obs*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*sin(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)) + 500*Y_obs*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*sin(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))^2)/250000 + (X_obs - cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*(sin(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*(Z_obs - 43/500) + Y_obs*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma)) + X_obs*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))))^2 + (Y_obs - cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(sin(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*(Z_obs - 43/500) + Y_obs*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma)) + X_obs*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))))^2)^(1/2);

            d2 = (((73*sin(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))/500 - Z_obs + sin(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*(cos(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(X_obs - (73*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))/500) - sin(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*((73*sin(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))/500 - Z_obs + 43/500) + cos(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(Y_obs - (73*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma)))/500)) + 43/500)^2 + ((73*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))/500 - X_obs + cos(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(cos(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(X_obs - (73*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))/500) - sin(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*((73*sin(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))/500 - Z_obs + 43/500) + cos(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(Y_obs - (73*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma)))/500)))^2 + ((73*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma)))/500 - Y_obs + cos(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(cos(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(X_obs - (73*cos(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))/500) - sin(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*((73*sin(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma)))/500 - Z_obs + 43/500) + cos(a10 + a15 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a11*(t_begin + p_tao*sigma)^4 + a12*(t_begin + p_tao*sigma)^3 + a13*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma) + a14*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma))*(Y_obs - (73*cos(a10 + a6*(t_begin + p_tao*sigma)^4 + a7*(t_begin + p_tao*sigma)^3 + a8*(t_begin + p_tao*sigma)^2 + a9*(t_begin + p_tao*sigma))*sin(a5 + a1*(t_begin + p_tao*sigma)^4 + a2*(t_begin + p_tao*sigma)^3 + a3*(t_begin + p_tao*sigma)^2 + a4*(t_begin + p_tao*sigma)))/500)))^2)^(1/2);
 
            diff1 = d1 - dtg.Threshold;
            diff2 = d2 - dtg.Threshold;
            f = 0;
            interval = 1/dtg.N;

            for tao = interval:interval:1

                    %%% solution from torque sparse regression without noise
                    e =  abs((4318789012341795*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/9007199254740992 + (943383695680083*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/1125899906842624 - (3954474643941301*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/576460752303423488 + (183003590714215*cos(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/36028797018963968 - (8970674701475047*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/2305843009213693952 + (858455264267093*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/144115188075855872 + (3209595118244909*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/288230376151711744 - (5628802016702833*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/288230376151711744 - (4707113331363697*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))^2*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1152921504606846976 + (1876048283815069*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/288230376151711744 - (4652447558222763*cos(2*a10 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))^2*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/576460752303423488) + abs((822897723913137*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^3)/1152921504606846976 + (79476643920073*cos(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/36028797018963968 - (5289765272146859*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/1152921504606846976 + (4793697736359671*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/1152921504606846976 - (4287075204197807*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/576460752303423488 + (1632836610114493*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/1152921504606846976 + (7265751674285011*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))^2)/576460752303423488 + (7252498841589555*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao))^2)/2305843009213693952 + (4679938971500113*cos(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/1152921504606846976 - (2889301995050081*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1152921504606846976 - (2403178407240127*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/576460752303423488) + abs((4110274331178377*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/9007199254740992 + (4990997193443041*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1152921504606846976 + (2517888332340985*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/576460752303423488 + (3784776127678219*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1152921504606846976 + (5741468388436775*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1152921504606846976 + (4555999909754877*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao))^2)/2305843009213693952 + (1681902067046799*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))^2*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/288230376151711744 + (8543932335759869*cos(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/4611686018427387904 - (2044694759205197*cos(2*a15 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao))^2)/1152921504606846976);
                  
                    %%%% analytical solution
%                    e = abs((a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao))*((34969*a8)/4000000 + (5140329351327242409843*a13)/576460752303423488000000 + (183447*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao)))/400000 + (104907*a6*(t_begin + sigma*tao)^2)/2000000 + (15420988053981727229529*a11*(t_begin + sigma*tao)^2)/288230376151711744000000 + (34969*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2)/16000000 + (13651*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao)))/2000000 + (13651*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2)/4000000 + (13651*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))^2)/2000000 + (104907*a7*(t_begin + sigma*tao))/4000000 + (15420988053981727229529*a12*(t_begin + sigma*tao))/576460752303423488000000 + (13651*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2)/4000000)) + abs((a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*((728679042008825725203*a8)/22517998136852480000000 + (34969*a13)/4000000 + (183447*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao)))/400000 + (214839*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/250000 + (2186037126026477175609*a6*(t_begin + sigma*tao)^2)/11258999068426240000000 + (104907*a11*(t_begin + sigma*tao)^2)/2000000 + (34969*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2)/16000000 + (13651*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao)))/1000000 + (13651*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao)))/2000000 - (13651*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao))^2)/2000000 + (2186037126026477175609*a7*(t_begin + sigma*tao))/22517998136852480000000 + (104907*a12*(t_begin + sigma*tao))/4000000 + (13651*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2)/2000000 + (58619*sin(2*a10 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2)/10000000 - (13651*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1000000)) + abs((a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*((5998407394089546817249*a3)/360287970189639680000000 + (17995222182268640451747*a1*(t_begin + sigma*tao)^2)/180143985094819840000000 + (13651*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao)))/2000000 + (13651*cos(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao)))/2000000 + (58619*cos(2*a10 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao)))/10000000 + (17995222182268640451747*a2*(t_begin + sigma*tao))/360287970189639680000000 + (34969*cos(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao)))/16000000 - (58619*sin(2*a10 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(2*a9 + 8*a6*(t_begin + sigma*tao)^3 + 6*a7*(t_begin + sigma*tao)^2 + 4*a8*(t_begin + sigma*tao)))/10000000 - (13651*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/2000000 - (34969*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(2*a9 + 2*a14 + 8*a6*(t_begin + sigma*tao)^3 + 6*a7*(t_begin + sigma*tao)^2 + 8*a11*(t_begin + sigma*tao)^3 + 6*a12*(t_begin + sigma*tao)^2 + 4*a8*(t_begin + sigma*tao) + 4*a13*(t_begin + sigma*tao)))/16000000 - (13651*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(2*a9 + a14 + 8*a6*(t_begin + sigma*tao)^3 + 6*a7*(t_begin + sigma*tao)^2 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 4*a8*(t_begin + sigma*tao) + 2*a13*(t_begin + sigma*tao)))/2000000));
                    
                    %%%  5% measurement noise, perfect fitting accuracy
%                    e =  abs((4086541802294015*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/9007199254740992 + (6343581644217701*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/2305843009213693952 + (3900598622031023*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1152921504606846976 + (6008842414570149*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/4611686018427387904 + (1706075948694643*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/576460752303423488 + (2445294629803415*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/288230376151711744 - (6877361242420579*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/2305843009213693952 + (2765040115283555*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1152921504606846976 + (6129380357876795*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/1152921504606846976) + abs((8593452656254539*cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/18014398509481984 + (7649089157439525*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/9007199254740992 + (1285455596168927*(2*a8 + 12*a6*(t_begin + sigma*tao)^2 + 6*a7*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/72057594037927936 + (2242541854003255*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/288230376151711744 - (6501762474649761*sin(2*a10 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/2305843009213693952 + (7537293511657539*cos(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/1152921504606846976 + (8513679675478985*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/2305843009213693952 + (5701611892022517*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/2305843009213693952 + (4889309516736717*cos(2*a10 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao))*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/576460752303423488 + (7038597314839847*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao))^2)/1152921504606846976 - (1631966154378515*cos(2*a10 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))^2*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/288230376151711744) + abs((2903102465460225*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/2305843009213693952 - (7168681448204637*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/2305843009213693952 - (4427022581034509*cos(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/576460752303423488 - (7965488558468521*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/2305843009213693952 + (4839537895382839*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/288230376151711744 - (1291478169878617*cos(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a3 + 12*a1*(t_begin + sigma*tao)^2 + 6*a2*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/144115188075855872 - (2624233811925921*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(2*a13 + 12*a11*(t_begin + sigma*tao)^2 + 6*a12*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao)))/2305843009213693952 + (6807194439650207*cos(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/576460752303423488 + (4030051530872041*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))^2*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao)))/288230376151711744 + (419208744258513*sin(2*a10 + 2*a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + 2*a11*(t_begin + sigma*tao)^4 + 2*a12*(t_begin + sigma*tao)^3 + 2*a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + 2*a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/72057594037927936 + (396505558104983*sin(a15 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/72057594037927936 + (1341957396805947*sin(2*a10 + a15 + 2*a6*(t_begin + sigma*tao)^4 + 2*a7*(t_begin + sigma*tao)^3 + 2*a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + 2*a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(a4 + 4*a1*(t_begin + sigma*tao)^3 + 3*a2*(t_begin + sigma*tao)^2 + 2*a3*(t_begin + sigma*tao))*(a9 + 4*a6*(t_begin + sigma*tao)^3 + 3*a7*(t_begin + sigma*tao)^2 + 2*a8*(t_begin + sigma*tao))*(a14 + 4*a11*(t_begin + sigma*tao)^3 + 3*a12*(t_begin + sigma*tao)^2 + 2*a13*(t_begin + sigma*tao)))/288230376151711744);

                    

                    f = f + e;
          
            end

            E_at = f * (interval*sigma) + dtg.Alpha * (max(0, -diff1)^2 + max(0, -diff2)^2);

            
        end
        
        
        
        
        function [c, c_eq] = myConstraints(dtg, a)

            t = dtg.T_end;
    
            c = [];
            c_eq = [t^4 * a(1) + t^3 * a(2) + t^2 * a(3) + t * dtg.d_theta_0(1) + dtg.theta_0(1) - dtg.theta_end(1);
                    t^4 * a(4) + t^3 * a(5) + t^2 * a(6) + t * dtg.d_theta_0(2) + dtg.theta_0(2) - dtg.theta_end(2);
                    t^4 * a(7) + t^3 * a(8) + t^2 * a(9) + t * dtg.d_theta_0(3) + dtg.theta_0(3) - dtg.theta_end(3);];


         end
        
         function [para] = gradientDescent(dtg)
             
                ObjFcn = @(a) dtg.objectFunction(a);

                ConsFcn = @(a) dtg.myConstraints(a);
                X = [];
                FVAL = [];
                for i = 0:5
                    options = optimset('display','final','TolX',1e-11,'MaxIter',100000,'MaxFunEvals',100000);
                    [x,fval] = fmincon(ObjFcn,rand(9,1),[],[],[],[],[],[],ConsFcn,options);
                    X = [X; x'];
                    FVAL = [FVAL; fval];
                end
                fval = min(FVAL);
                selected = FVAL == fval;
                para = X' * selected;
        
         end
         
         
         

            function [theta] = joint_space(dtg,t,a)

            a5 = dtg.theta_0(1);
            a10 = dtg.theta_0(2);
            a15 = dtg.theta_0(3);

            a4 = dtg.d_theta_0(1);
            a9 = dtg.d_theta_0(2);
            a14 = dtg.d_theta_0(3);
            
            theta = zeros(3, 1);
            theta(1) = a(1)*t^4 + a(2)*t^3 + a(3)*t^2 + a4*t +a5;
            theta(2) = a(4)*t^4 + a(5)*t^3 + a(6)*t^2 + a9*t +a10;
            theta(3) = a(7)*t^4 + a(8)*t^3 + a(9)*t^2 + a14*t +a15;

            end

            
            
            function [M] = geForce(dtg,t, a)


            a5 = dtg.theta_0(1);
            a10 = dtg.theta_0(2);
            a15 = dtg.theta_0(3);


            a4 = dtg.d_theta_0(1);
            a9 = dtg.d_theta_0(2);
            a14 = dtg.d_theta_0(3);

            a1 = a(1);
            a2 = a(2);
            a3 = a(3);
            a6 = a(4);
            a7 = a(5);
            a8 = a(6) ;
            a11 = a(7);
            a12 = a(8);
            a13 = a(9) ;

            M = zeros(3,1);

            M(1) =  (76779614644346195*a3)/4611686018427387904 + (3934632864847017*cos(2*a10 + a15 + 2*a9*t + a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + a11*t^4 + a12*t^3 + a13*t^2)*(12*a1*t^2 + 6*a2*t + 2*a3))/576460752303423488 + (230338843933038585*a2*t)/4611686018427387904 + (3934632864847017*cos(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(12*a1*t^2 + 6*a2*t + 2*a3))/576460752303423488 + (1259891002956151*cos(2*a10 + 2*a15 + 2*a9*t + 2*a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + 2*a11*t^4 + 2*a12*t^3 + 2*a13*t^2)*(12*a1*t^2 + 6*a2*t + 2*a3))/576460752303423488 + (230338843933038585*a1*t^2)/2305843009213693952 + (27033242271419503*cos(2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + 2*a9*t + 2*a10)*(12*a1*t^2 + 6*a2*t + 2*a3))/4611686018427387904 - (27033242271419503*sin(2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + 2*a9*t + 2*a10)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)*(4*a6*t^3 + 3*a7*t^2 + 2*a8*t + a9))/2305843009213693952 - (3934632864847017*sin(2*a10 + a15 + 2*a9*t + a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + a11*t^4 + a12*t^3 + a13*t^2)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)*(4*a6*t^3 + 3*a7*t^2 + 2*a8*t + a9))/288230376151711744 - (3934632864847017*sin(2*a10 + a15 + 2*a9*t + a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + a11*t^4 + a12*t^3 + a13*t^2)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)*(4*a11*t^3 + 3*a12*t^2 + 2*a13*t + a14))/576460752303423488 - (3934632864847017*sin(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)*(4*a11*t^3 + 3*a12*t^2 + 2*a13*t + a14))/576460752303423488 - (1259891002956151*sin(2*a10 + 2*a15 + 2*a9*t + 2*a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + 2*a11*t^4 + 2*a12*t^3 + 2*a13*t^2)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)*(4*a6*t^3 + 3*a7*t^2 + 2*a8*t + a9))/288230376151711744 - (1259891002956151*sin(2*a10 + 2*a15 + 2*a9*t + 2*a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + 2*a11*t^4 + 2*a12*t^3 + 2*a13*t^2)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)*(4*a11*t^3 + 3*a12*t^2 + 2*a13*t + a14))/288230376151711744;

            M(2) = (37308366950851875*a8)/1152921504606846976 + (1259891002956151*a13)/144115188075855872 + (48377442017232500631*cos(a6*t^4 + a7*t^3 + a8*t^2 + a9*t + a10))/56294995342131200000 + (4130859204211177*cos(a10 + a15 + a9*t + a14*t + a6*t^4 + a7*t^3 + a8*t^2 + a11*t^4 + a12*t^3 + a13*t^2))/9007199254740992 + (1259891002956151*sin(2*a10 + 2*a15 + 2*a9*t + 2*a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + 2*a11*t^4 + 2*a12*t^3 + 2*a13*t^2)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)^2)/576460752303423488 + (111925100852555625*a7*t)/1152921504606846976 + (3779673008868453*a12*t)/144115188075855872 + (3934632864847017*cos(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(12*a6*t^2 + 6*a7*t + 2*a8))/288230376151711744 + (3934632864847017*cos(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(12*a11*t^2 + 6*a12*t + 2*a13))/576460752303423488 + (27033242271419503*sin(2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + 2*a9*t + 2*a10)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)^2)/4611686018427387904 + (111925100852555625*a6*t^2)/576460752303423488 + (3779673008868453*a11*t^2)/72057594037927936 + (3934632864847017*sin(2*a10 + a15 + 2*a9*t + a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + a11*t^4 + a12*t^3 + a13*t^2)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)^2)/576460752303423488 - (3934632864847017*sin(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(4*a11*t^3 + 3*a12*t^2 + 2*a13*t + a14)^2)/576460752303423488 - (3934632864847017*sin(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(4*a6*t^3 + 3*a7*t^2 + 2*a8*t + a9)*(4*a11*t^3 + 3*a12*t^2 + 2*a13*t + a14))/288230376151711744;

            M(3) =  (1259891002956151*a8)/144115188075855872 + (328981078484943515*a13)/36893488147419103232 + (4130859204211177*cos(a10 + a15 + a9*t + a14*t + a6*t^4 + a7*t^3 + a8*t^2 + a11*t^4 + a12*t^3 + a13*t^2))/9007199254740992 + (1259891002956151*sin(2*a10 + 2*a15 + 2*a9*t + 2*a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + 2*a11*t^4 + 2*a12*t^3 + 2*a13*t^2)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)^2)/576460752303423488 + (3779673008868453*a7*t)/144115188075855872 + (986943235454830545*a12*t)/36893488147419103232 + (3934632864847017*cos(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(12*a6*t^2 + 6*a7*t + 2*a8))/576460752303423488 + (3779673008868453*a6*t^2)/72057594037927936 + (986943235454830545*a11*t^2)/18446744073709551616 + (3934632864847017*sin(2*a10 + a15 + 2*a9*t + a14*t + 2*a6*t^4 + 2*a7*t^3 + 2*a8*t^2 + a11*t^4 + a12*t^3 + a13*t^2)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)^2)/1152921504606846976 + (3934632864847017*sin(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4)^2)/1152921504606846976 + (3934632864847017*sin(a11*t^4 + a12*t^3 + a13*t^2 + a14*t + a15)*(4*a6*t^3 + 3*a7*t^2 + 2*a8*t + a9)^2)/576460752303423488;


            end

            
            
            function [v] = geVelo(dtg,t, a)

            a4 = dtg.d_theta_0(1);
            a9 = dtg.d_theta_0(2);
            a14 = dtg.d_theta_0(3);


            a1 = a(1);
            a2 = a(2);
            a3 = a(3);
            a6 = a(4);
            a7 = a(5);
            a8 = a(6) ;
            a11 = a(7);
            a12 = a(8);
            a13 = a(9) ;


            v = zeros(3,1);
            v(1) = 4*a1*t^3 + 3*a2*t^2 + 2*a3*t + a4;
            v(2) = 4*a6*t^3 + 3*a7*t^2 + 2*a8*t + a9;
            v(3) = 4*a11*t^3 + 3*a12*t^2 + 2*a13*t + a14;

            end

            
            
            function [A] = geAcce(dtg,t, a)

            a1 = a(1);
            a2 = a(2);
            a3 = a(3);
            a6 = a(4);
            a7 = a(5);
            a8 = a(6) ;
            a11 = a(7);
            a12 = a(8);
            a13 = a(9) ;

            A = zeros(3,1);
            A(1) = 12*a1*t^2 + 6*a2*t + 2*a3;
            A(2) = 12*a6*t^2 + 6*a7*t + 2*a8;
            A(3) = 12*a11*t^2 + 6*a12*t + 2*a13;

            end
            
            
            
           function [e,e_sp] = energy(dtg,I_t,a)

            M = geForce(I_t,a,dtg.theta_0, dtg.d_theta_0);
            v = geVelo(I_t, a, dtg.theta_0, dtg.d_theta_0);

            e_sp = M.*v;
            e = sum(abs(e_sp));

           end
           
                

            function [result] = getResult(dtg)
                
                x = dtg.gradientDescent();
                %%point_number = 40, 0-40,41 points
                path = zeros(dtg.Point_number+1,3);
                CarPoints = zeros(dtg.Point_number+1, 3);
                path(1,:) =  dtg.joint_space(0, x);
                CarMat = dtg.lynxmotion.fkine(path(1,:)).T;
                CarPoints(1, :) = CarMat(1:3, 4, 1);
                VELO = zeros(dtg.Point_number+1,3);
                VELO(1,:) = dtg.geAcce(0, x);


                for i = 1:dtg.Point_number

                    path(i+1,:) =  dtg.joint_space((1/dtg.Point_number)*dtg.Sigma*i, x);
                    CarMat = dtg.lynxmotion.fkine(path(i+1,:)).T;
                    CarPoints(i+1, :) = CarMat(1:3, 4);
                    VELO(i+1,:) = dtg.geAcce(i*(1/dtg.Point_number)*dtg.Sigma, x);
                end
                
                result.a = x;
                dtg.parameter = x;
                result.Path = path;
                result.Carpoints = CarPoints;
                result.Velo = VELO;

     
                if strcmp(dtg.state,'initialization')
                    dtg.global_path = path;
                    dtg.global_carpoints = CarPoints;
                    dtg.global_velo = VELO;
                    d = dtg.distance();
                    dtg.global_distance = [d.d1 d.d2];
                    dtg.init_distance = [d.d1 d.d2];
                    dtg.init_path = path;
                    dtg.init_velo = VELO;
                    dtg.init_carpoints = CarPoints;
                    for n = 1:dtg.global_point_number
                        dtg.time = n;
                        dtg.P_tao = n/dtg.global_point_number;      
                        init_d = dtg.distance();
                        dtg.init_distance = [dtg.init_distance; init_d.d1 init_d.d2];
                    end
                    dtg.P_tao = 0;
                    dtg.time = 0;
                    
                    
                elseif strcmp(dtg.state, 'dtg')
                    
                    dtg.global_path = [dtg.global_path(1:dtg.time,:); path];
                    dtg.global_carpoints = [dtg.global_carpoints(1:dtg.time,:); CarPoints];
                    dtg.global_velo = [dtg.global_velo(1:dtg.time,:); VELO];
                    
                end

            end
            
            
            
           function [D] = distance(dtg)
               
                D.a = dtg.parameter;

                a1 = D.a(1);
                a2 = D.a(2);
                a3 = D.a(3);
                a6 = D.a(4);
                a7 = D.a(5);
                a8 = D.a(6);
                a11 = D.a(7);
                a12 = D.a(8);
                a13 = D.a(9);

                a5 = dtg.theta_0(1);
                a10 = dtg.theta_0(2);
                a15 = dtg.theta_0(3);


                a4 = dtg.d_theta_0(1);
                a9 = dtg.d_theta_0(2);
                a14 = dtg.d_theta_0(3);

                X_obs = dtg.Object(dtg.time+1,1);
                Y_obs = dtg.Object(dtg.time+1,2);
                Z_obs = dtg.Object(dtg.time+1,3);

                t_begin = dtg.T_begin;
                tao = dtg.P_tao;
                sigma = dtg.Sigma;



                D.d1 = ((cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))^2*(43*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)) - 500*Z_obs*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)) + 500*X_obs*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*sin(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)) + 500*Y_obs*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*sin(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))^2)/250000 + (X_obs - cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*(sin(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*(Z_obs - 43/500) + X_obs*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)) + Y_obs*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))))^2 + (Y_obs - cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(sin(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*(Z_obs - 43/500) + X_obs*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)) + Y_obs*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))))^2)^(1/2);

                D.d2 = (((73*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao)))/500 - Y_obs + cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(X_obs - (73*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/500) - sin(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*((73*sin(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/500 - Z_obs + 43/500) + cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(Y_obs - (73*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao)))/500)))^2 + ((73*sin(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/500 - Z_obs + sin(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*(cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(X_obs - (73*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/500) - sin(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*((73*sin(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/500 - Z_obs + 43/500) + cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(Y_obs - (73*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao)))/500)) + 43/500)^2 + ((73*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/500 - X_obs + cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(X_obs - (73*cos(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/500) - sin(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*((73*sin(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao)))/500 - Z_obs + 43/500) + cos(a10 + a15 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a11*(t_begin + sigma*tao)^4 + a12*(t_begin + sigma*tao)^3 + a13*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao) + a14*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao))*(Y_obs - (73*cos(a10 + a6*(t_begin + sigma*tao)^4 + a7*(t_begin + sigma*tao)^3 + a8*(t_begin + sigma*tao)^2 + a9*(t_begin + sigma*tao))*sin(a5 + a1*(t_begin + sigma*tao)^4 + a2*(t_begin + sigma*tao)^3 + a3*(t_begin + sigma*tao)^2 + a4*(t_begin + sigma*tao)))/500)))^2)^(1/2);

                D.diff1 = D.d1 - dtg.Threshold;
                D.diff2 = D.d2 - dtg.Threshold;
                
                
            end
            
            
            
            function dynamicTrajGen(dtg)
                
                time_avoidance = 0;
                dtg.state = 'dtg';
                for k = 1:dtg.global_point_number
                    dtg.time = k;    
                    dtg.P_tao = (1/dtg.Point_number)*(dtg.time-time_avoidance);
                    d = dtg.distance();
                    d.d1
                    d.d2
                    dtg.global_distance = [dtg.global_distance; d.d1 d.d2];

                    if d.d1<(dtg.Threshold)||d.d2<(dtg.Threshold)

                    time_avoidance = dtg.time;  


                    dtg.T_begin = 1/(dtg.global_point_number) * dtg.global_sigma * dtg.time;
                    dtg.Sigma = dtg.T_end - dtg.T_begin;

                    dtg.theta_0 = dtg.global_path(dtg.time+1,:);
                    dtg.d_theta_0 = dtg.global_velo(dtg.time+1,:);
                    dtg.Point_number = dtg.global_point_number-dtg.time;

                    dtg.P_tao = 0;

                    dtg.getResult();
 
                end

                end

            end
            
             
            
            function plot_distance(dtg)
                
                figure(1);
                T = linspace(0,200,201);
                thre = zeros(201,1)+0.08;
                hold on;
                plot(T,dtg.global_distance(:,1),'-r','LineWidth',1.5 );
                plot(T,dtg.global_distance(:,2),'-k','LineWidth',1.5 );
                plot(T,dtg.init_distance(:,1),'--m','LineWidth',1.5 );
                plot(T,dtg.init_distance(:,2),'--b','LineWidth',1.5 );
                plot(T,thre,'-g','LineWidth',1.5 );
                xlabel('Time(1/100s)');
                ylabel('Distance to obstacle(m)');
                legend('d1 - with obstacle avoidance', 'd2 - with obstacle avoidance', 'd1 - without obstacle avoidance', 'd2 - without obstacle avoidance', 'threshold')
                
            end

            
            function plot_robot_trajectory(dtg,P)
                
                figure(2);
                hold on;
                text( 0.15+0.02,    0.4+0.02,   0.22+0.02, 'ObstacleStart');
                text( 0.18-0.01,    -0.35-0.08,   0.28+0.02, 'ObstacleEnd');

                text( 0.0+0.02,    0.3028-0.03,   0.1426-0.02, 'RobotStart');
                text( 0.0-0.01,    -0.228-0.08,   0.2637+0.03, 'RobotEnd');

                for i = 1:dtg.global_point_number+1
                    dtg.lynxmotion.plot(dtg.global_path(i,:),'jointdiam',1,'trail','b-')
                    plot3(P(1:i,1), P(1:i,2), P(1:i,3),'-k','LineWidth',3 );
                    plot3(dtg.init_carpoints(1:i,1), dtg.init_carpoints(1:i,2), dtg.init_carpoints(1:i,3),'-r','LineWidth',3);   
                    plot3(dtg.global_carpoints(1:i,1), dtg.global_carpoints(1:i,2), dtg.global_carpoints(1:i,3),'-g','LineWidth',3);   
                end

                xlabel('X(m)')
                ylabel('Y(m)')
                zlabel('Z(m)') 
               
            end
            

            
         
         
         
         
         
    end
    
    
    
    
end