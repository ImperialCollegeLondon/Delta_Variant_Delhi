data 
{
  int <lower=1> M;               // number of regions
  int <lower=1> N0;              // number of days in init period
  int<lower=1> N[M];             // number of days in entire period
  int<lower=1> N2;               // max number of days across regions
  int deaths[N2, M];             // observed deaths 
  matrix[N2, M] f;               // infection-to-death distribution
  int ll_len;                    // length of indices for likelihood
  int ll_idxs[ll_len];           // indices for likelihood
  real pop[M];                   // population count
  int W;                         // maximum number of weeks
  int week_index[M,N2];
  real SI[N2];                   // generation time distribution
  real WI[N2];                   // waning immunity distribution
  real AR_SD_MEAN;
  int T2;                        // starting index of the second strain
  int phylo_N_len; 
  int phylo_N[phylo_N_len]; 
  int phylo_PSamples[phylo_N_len];
  int phylo_NSamples[phylo_N_len];
  int sero_dates_len;
  int sero_dates[sero_dates_len];     // dates on which the serological surveys were done
  int sero_tested[sero_dates_len];  // number of individuals tested
  int sero_pos[sero_dates_len]; // number of positive samples
  matrix[N2, M] PCR_pos_prob;     // PCR positivity distribution
  matrix[N2, M] seroconv_cdf;     // cumululative seroconverted distribution
  matrix[N2, M] serorev_surv;  // cumululative seroreverted distribution
  real <lower=0,upper=1> UR; // underreporting factor
  real ifrMean; // ifr mean to use
  real ifrSD; // ifr sd to use
  real RdiffPar1; // parameter 1 of gamma distribution for R_difference
  real RdiffPar2; // parameter 2 of gamma distribution for R_difference
  real crossPar1; // parameter 2 of beta distribution for cross
  real crossPar2; // parameter 2 of beta distribution for cross
  int max_delay; //  
}

transformed data {
  int i_delay_mask[N2];
  vector[N2] SI_rev; // SI in reverse order
  vector[N2] WI_rev; // waning immunity in revered order
  matrix[N2, M] seroconv_cdf_rev; // cumululative seroconverted distribution in reverse order
  matrix[N2, M] serorev_surv_rev; // cumululative seroreverted distribution in reverse order
  matrix[N2, M] serorev_conv_rev; // cumululative seroreverted and seroconverted distribution in reverse order
  matrix[N2, M] PCR_pos_prob_rev; // PCR positivity distribution in reverse order
  matrix[N2, M] f_rev; // infection-to-death distribution in reversed order
  for(i in 1:N2){
    SI_rev[i] = SI[N2-i+1];
    WI_rev[i] = WI[N2-i+1];
    seroconv_cdf_rev[i,:] = seroconv_cdf[N2-i+1,:];
    serorev_surv_rev[i,:] = serorev_surv[N2-i+1,:];
    PCR_pos_prob_rev[i,:] = PCR_pos_prob[N2-i+1,:];
    f_rev[i,:] = f[N2-i+1,:];
    serorev_conv_rev[i,:] = seroconv_cdf_rev[i,:] .* serorev_surv_rev[i,:];
    if (i<max_delay){
      i_delay_mask[i] = i;
    }
    else{
      i_delay_mask[i] = max_delay;
    }
  }
  
}

parameters 
{
  real<lower=0> R_difference; // relative transmissibility of strain 2 compared to strain 1 
  real<lower=1> y_v1[M]; // infections in intial period strain 1
  real<lower=1> y_v2[M]; // infections in intial period strain 2
  real<lower=0> phi;
  real<lower=0,upper=1> cross; // cross-immunity
  real<lower=0> tau; // NegBin dispersion
  matrix[W+1,M] weekly_effect;
  real<lower=0, upper=1> weekly_rho;
  real<lower=0, upper=1> weekly_rho1;
  real<lower=0> weekly_sd;
  real <lower=0> ifr1[M];
  real <lower=0> ifr2[M];
  real sero_effect_mu1;
  real<lower=0> sero_delta1;
  real sero_effect1[(sero_dates_len-2)];
}


transformed parameters 
{
  matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
  matrix[N2,M] prediction = rep_matrix(0,N2,M);
  matrix[N2,M] E_deaths  = rep_matrix(0,N2,M);
  matrix[N2,M] immune = rep_matrix(0,N2,M);
  matrix[N2,M] prediction_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] prediction_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] E_deaths_v1  = rep_matrix(0,N2,M);
  matrix[N2,M] E_deaths_v2  = rep_matrix(0,N2,M);
  matrix[N2,M] Rt_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] Rt_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] Rt_adj_immune_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] Rt_adj_immune_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] cumm_sum_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] cumm_sum_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] immune_v1 = rep_matrix(0,N2,M);
  matrix[N2,M] immune_v2 = rep_matrix(0,N2,M);
  matrix[N2,M] alpha_sus1 = rep_matrix(0,N2,M);
  matrix[N2,M] alpha_sus2 = rep_matrix(0,N2,M);
  matrix[N2,M] ar1 = rep_matrix(0,N2,M);
  matrix[N2,M] ar2 = rep_matrix(0,N2,M);
  matrix[N2,M] arboth = rep_matrix(0,N2,M);
  matrix[N2,M] ar = rep_matrix(0,N2,M);
  matrix[N2,M] n1 = rep_matrix(0,N2,M);
  matrix[N2,M] n2 = rep_matrix(0,N2,M);
  matrix[N2,M] seropos_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] seroconv_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] seroconv_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] pcr_pos_v1 = rep_matrix(0, N2, M);
  matrix[N2,M] pcr_pos_v2 = rep_matrix(0, N2, M);
  matrix[N2,M] E_fraction = rep_matrix(0,N2,M);
  real <lower=0> RR[M]; // relative risk
  matrix[N2,M] seroconv_v1_only_fraction = rep_matrix(0, N2, M);
  matrix[N2,M] seroconv_v2_only_fraction = rep_matrix(0, N2, M);
  matrix[N2,M] seroconv_both_fraction = rep_matrix(0, N2, M);
  matrix[N2,M] seroconv_fraction = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_v1_only_fraction = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_v2_only_fraction = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_both_fraction = rep_matrix(0, N2, M);
  matrix[N2,M] seropos_fraction = rep_matrix(0, N2, M);  

  for (m in 1:M)
  {
    for (i in 2:N0)
    {
      cumm_sum_v1[i,m] = cumm_sum_v1[i-1,m] + y_v1[m];
    }
    
    for ( i in ( T2 + 1 ):( T2 + 1 ) )
    {
      cumm_sum_v2[i,m] = cumm_sum_v2[i-1,m] + y_v2[m];
    }

    prediction_v1[1:N0,m] = rep_vector(y_v1[m],N0); 
    prediction_v2[T2:(T2+N0-1),m] = rep_vector(y_v2[m],N0); // strain2 
    
    // reproduction numbers
    Rt_v1[,m] = 3.28 * 2 * inv_logit( - weekly_effect[week_index[m],m] );
    Rt_v2[T2:N2,m] = Rt_v1[T2:N2,m] * R_difference; 
    
    // adjusted reproduction number during initial period
    Rt_adj_immune_v1[1:N0,m] = Rt_v1[1:N0,m]; 
    
    for (i in (N0+1):N2) 
    {
      
      // strain 1 //
      real convolution_v1 = dot_product(tail(sub_col(prediction_v1, 1, m, i-1), i_delay_mask[i]-1), tail(SI_rev, i_delay_mask[i]-1));
      immune_v1[i,m] = dot_product(sub_col(prediction_v1, 1, m, i-1), tail(WI_rev, i-1));

      // strain 2 //
      if ( i > T2 ) 
      {
        real convolution_v2 = 0; 
        convolution_v2 = dot_product(tail(sub_col(prediction_v2, 1, m, i-1), i_delay_mask[i]-1), tail(SI_rev, i_delay_mask[i]-1));
        immune_v2[i,m] = dot_product(sub_col(prediction_v2, 1, m, i-1 ), tail(WI_rev, i-1 ));
        alpha_sus2[i,m] = ( 1 - cross ) * immune_v2[i,m] / ( pop[m] - cross * immune_v2[i,m] );
        n2[i,m] = immune_v2[i,m] + cross * ( immune_v1[i,m] * ( 1 - alpha_sus2[i,m] ) );
        Rt_adj_immune_v2[i,m] = ( 1 - n2[i,m] / pop[m]) * Rt_v2[i,m];
        prediction_v2[i, m] = ( pop[m] - n2[i,m] ) * ( 1 - exp( -Rt_v2[i,m] * convolution_v2 / pop[m] ) );
        cumm_sum_v2[i,m]  = cumm_sum_v2[i-1,m] +  prediction_v2[i,m];
      }
      alpha_sus1[i,m] = (1 - cross) * immune_v1[i,m] / ( pop[m] - cross * immune_v1[i,m] );
      n1[i,m] = immune_v1[i,m] + cross * ( immune_v2[i,m] * ( 1 - alpha_sus1[i,m] ) );
      Rt_adj_immune_v1[i,m] = ( 1 - n1[i,m] / pop[m]) * Rt_v1[i,m];
      prediction_v1[i, m] = ( pop[m] - n1[i,m] ) * ( 1 - exp( -Rt_v1[i,m] * convolution_v1 / pop[m] ));
      cumm_sum_v1[i,m] = cumm_sum_v1[i-1,m] +  prediction_v1[i,m];
      
      cumm_sum[i,m] = cumm_sum_v1[i,m] + cumm_sum_v2[i,m];
      prediction[i, m] = prediction_v1[i, m] + prediction_v2[i, m];
      immune[i, m] = immune_v1[i, m] + immune_v2[i, m];
    }
    
    E_deaths_v1[1, m]= 1e-15 * prediction_v1[1,m];
    E_deaths_v2[1, m]= 1e-15 * prediction_v2[1,m];
    E_deaths[1, m]= 1e-15 * (prediction_v1[1,m] + prediction_v2[1,m]);
    
    for (i in 2:N2)
    {
      // strain 1 //
      // number of positive PCR 
      pcr_pos_v1[i,m] = dot_product(sub_col(prediction_v1, 1, m, i-1), tail(PCR_pos_prob_rev[:,m], i-1));
      
      // number of seroconverted and reverted 
      seropos_v1[i,m] = 0.486 * dot_product(sub_col(prediction_v1, 1, m, i-1), tail(serorev_conv_rev[:,m], i-1)) +
      0.514 * dot_product(sub_col(prediction_v1, 1, m, i-1), tail(seroconv_cdf_rev[:,m], i-1));
      seroconv_v1[i,m] = dot_product(sub_col(prediction_v1, 1, m, i-1), tail(seroconv_cdf_rev[:,m], i-1));
      
      // new deaths
      E_deaths_v1[i,m] = dot_product(tail(sub_col(prediction_v1, 1, m, i-1), i_delay_mask[i]-1), tail(f_rev[:,m], i_delay_mask[i]-1));
      E_deaths_v1[i,m] *= ifr1[m];
      
      
      // strain 2 //
      if (i > T2)  
      {
        // number of positive PCR 
        pcr_pos_v2[i,m] = dot_product(sub_col(prediction_v2, 1, m, i-1), tail(PCR_pos_prob_rev[:,m], i-1));
        
        // number of seroconverted and reverted 
        seropos_v2[i,m] = 0.486 * dot_product(sub_col(prediction_v2, 1, m, i-1), tail(serorev_conv_rev[:,m], i-1)) +
        0.514 * dot_product(sub_col(prediction_v2, 1, m, i-1), tail(seroconv_cdf_rev[:,m], i-1));
        seroconv_v2[i,m] = dot_product(sub_col(prediction_v2, 1, m, i-1), tail(seroconv_cdf_rev[:,m], i-1));
        
        // new deaths
        E_deaths_v2[i,m] = dot_product(tail(sub_col(prediction_v2, 1, m, i-1), i_delay_mask[i]-1), tail(f_rev[:,m], i_delay_mask[i]-1));
        E_deaths_v2[i,m] *= ifr2[m];
      }

      seroconv_v1_only_fraction[i,m] = (seroconv_v1[i,m] / pop[m]) * (1 - seroconv_v2[i,m] / pop[m]);
      seroconv_v2_only_fraction[i,m] = (seroconv_v2[i,m] / pop[m]) * (1 - seroconv_v1[i,m] / pop[m]);
      seroconv_both_fraction[i,m] = (seroconv_v1[i,m] / pop[m]) * (seroconv_v2[i,m] / pop[m]);
      seroconv_fraction[i,m] = seroconv_v1_only_fraction[i,m] + seroconv_v2_only_fraction[i,m] + seroconv_both_fraction[i,m];

      seropos_v1_only_fraction[i,m] = (seropos_v1[i,m] / pop[m]) * (1 - seropos_v2[i,m] / pop[m]);
      seropos_v2_only_fraction[i,m] = (seropos_v2[i,m] / pop[m]) * (1 - seropos_v1[i,m] / pop[m]);
      seropos_both_fraction[i,m] = (seropos_v1[i,m] / pop[m]) * (seropos_v2[i,m] / pop[m]);
      seropos_fraction[i,m] = seropos_v1_only_fraction[i,m] + seropos_v2_only_fraction[i,m] + seropos_both_fraction[i,m];
    }
    
    // sum of new cases, immunes and new deaths over the strains
    E_deaths[:,m] = E_deaths_v1[:,m] + E_deaths_v2[:,m];
    RR[m] = ifr2[m] / ifr1[m];
    
    // infection ratio
    E_fraction[1,m] = 0.0;
    E_fraction[2:N2,m] = pcr_pos_v2[2:N2,m]./(pcr_pos_v1[2:N2,m] + pcr_pos_v2[2:N2,m]);	
  }
}


model 
{
  ifr1 ~ normal(ifrMean,ifrSD);
  ifr2 ~ normal(ifrMean,ifrSD);
  cross ~ beta(crossPar1,crossPar2);
  R_difference ~ gamma(RdiffPar1,RdiffPar2);

  tau ~ exponential(0.03);
  phi ~ normal(0,5);

  sero_effect_mu1 ~ normal(0, 1);
  sero_delta1 ~ normal(0,1);
  sero_effect1 ~ normal(0,sero_delta1);

  weekly_sd ~ normal(0, AR_SD_MEAN);
  weekly_rho ~ normal(0.8, 0.05);
  weekly_effect[1, :] ~ normal(0, 0.01);
  weekly_effect[2, :] ~ normal(0,weekly_sd *sqrt(1-pow(weekly_rho, 2)));

  for (m in 1:M)
  {
    y_v1[m] ~ exponential(1/tau); 
    y_v2[m] ~ normal(0,20); 
    weekly_effect[3:(W+1), m] ~ normal( weekly_effect[2:W,m] * weekly_rho,weekly_sd * sqrt(1-pow(weekly_rho,2)) );
  }

  for (m in 1:M)
  {
    deaths[ll_idxs, m] ~ neg_binomial_2(E_deaths[ll_idxs, m] * UR, phi);
    for ( i in 1:phylo_N_len )
    {
      phylo_PSamples[i] ~ binomial(phylo_NSamples[i]+phylo_PSamples[i], E_fraction[phylo_N[i],m]);
    }
    for ( i in 1:(sero_dates_len-2) ) // sero survey 1
    { 
       sero_pos[i] ~ binomial(sero_tested[i], seropos_fraction[sero_dates[i],m] * 2 * inv_logit( sero_effect_mu1 + sero_effect1[i] ));
    }
    for ( i in (sero_dates_len-1):sero_dates_len ) // sero survey 2
    {
       sero_pos[i] ~ binomial(sero_tested[i], seroconv_fraction[sero_dates[i],m] ); // no sero reversion, no RE (not identiable)
    }
  }
}



generated quantities
{
  real cross_prior = beta_rng(crossPar1,crossPar2);
  real R_difference_prior = gamma_rng(RdiffPar1,RdiffPar2);
  real ifr1_prior = normal_rng(ifrMean,ifrSD);
  real ifr2_prior = normal_rng(ifrMean,ifrSD);
  real sero_effect_mu1_prior = normal_rng(0,1);
  real sero_delta1_prior = normal_rng(0,1);
  real sero_effect1_prior = normal_rng(0,sero_delta1);
}
