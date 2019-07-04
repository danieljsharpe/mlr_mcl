// original regularise (before parallelisation)
MLR_MCL::Csr_mtx MLR_MCL::regularise(const Csr_mtx &T_csr, const Csr_mtx_struct &TG_csr) {

    Csr_mtx T_csr_new;
    T_csr_new.first.reserve(n_row_occ*T_csr.second.size());
    T_csr_new.second.resize(T_csr.second.size());
    for (int i=0;i<T_csr.second.size();i++) { T_csr_new.second[i]=0; }
    signed int lo=0, hi=T_csr.second[0]; // low/high indices of elements of current row in transition matrix
    cout << "in regularise()" << endl;
    int rn=0;
    const vector<pair<double,int>> *const vp_T_er = &T_csr.first; // ptr to vec of elems & col indices of transn mtx
    const vector<double> *const vp_TG_e = &TG_csr.T_e;
    const vector<int> *const vp_TG_r = &TG_csr.T_c, *const vp_TG_cl = &TG_csr.T_rl;
    for (int i=0;i<n_row_occ;i++) {
        while (lo==hi) {
            lo = hi; hi = T_csr.second[rn+1]; rn++; }
        // low/high indices of elements of current col in regularisation matrix
        signed int lo_col=0, hi_col=(*vp_TG_cl)[0];
        int cn=0;
        for (int j=0;j<T_csr.second.size();j++) {
            if (lo_col==hi_col) {
                cout << "Error: column " << j+1 << " of regularisation matrix is empty" << endl;
                exit(EXIT_FAILURE); }
            int curr_row=lo, curr_col=lo_col;
            double val=0.;
            // multiplication of two sparse vectors with O(m+n) complexity
            while ((curr_row < hi) && (curr_col < hi_col)) {
                if ((*vp_T_er)[curr_row].second < (*vp_TG_r)[curr_col]) { curr_row++;
                } else if ((*vp_T_er)[curr_row].second > (*vp_TG_r)[curr_col]) { curr_col++;
                } else {
                    val += (*vp_T_er)[curr_row].first*(*vp_TG_e)[curr_col];
                    curr_row++; curr_col++;
                }
            }
            if (val>eps) { T_csr_new.first.emplace_back(make_pair(val,cn)); T_csr_new.second[rn]++; }
            lo_col = hi_col; hi_col = (*vp_TG_cl)[cn+1]; cn++;
        }
        lo = hi; hi = T_csr.second[rn+1]; rn++;
    }
    for (int i=1;i<T_csr_new.second.size();i++) {
        T_csr_new.second[i] += T_csr_new.second[i-1]; }
    T_csr_new.first.shrink_to_fit();
    cout << "leaving regularise()" << endl;
    return T_csr_new;
}

// alternative version of regularise function (note: currently does not work, has minor bug)
MLR_MCL::Csr_mtx MLR_MCL::regularise_alt(const Csr_mtx &T_csr, const Csr_mtx_struct &TG_csr) {

    cout << "in regularise_alt()..." << endl;
    vector<pair<double,int>>::const_iterator it_vec;
    const vector<double> *const vp_TG_e = &TG_csr.T_e;
    const vector<int> *const vp_TG_r = &TG_csr.T_c, *const vp_TG_cl = &TG_csr.T_rl;
    vector<vector<pair<double,int>>> T_ec_new(T_csr.second.size());
    vector<int> T_rl(T_csr.second.size(),0);
    int T_csr_nelems = T_csr.first.size();
    int m=0; // keep track of elems of TG_csr being scanned over
    vector<double> curr_colvec(T_csr.second.size(),0.); // column vector of current iteration
    for (int i=0;i<T_csr.second.size();i++) {
        int cn_lim=(*vp_TG_cl)[i];
        fill(curr_colvec.begin(),curr_colvec.end(),0.);
        do {
            curr_colvec[(*vp_TG_r)[m]] = (*vp_TG_e)[m]; m++;
        } while (m<cn_lim);
        int rn=0, k=0, rn_lim = T_csr.second[rn];
        double val=0.;
        if (k==rn_lim) while (k==rn_lim) { rn++; rn_lim = T_csr.second[rn]; }
        for (it_vec=T_csr.first.begin();it_vec!=T_csr.first.end();it_vec++) {
            val += it_vec->first*curr_colvec[it_vec->second];
            k++;
            if ((k==rn_lim) && (k<T_csr_nelems)) {
                if (val>eps) { T_ec_new[rn].emplace_back(make_pair(val,i)); T_rl[rn]++; }
                val=0.;
                while (k==rn_lim) { rn++; rn_lim = T_csr.second[rn]; };
            }
        }
    }
    vector<pair<double,int>> T_ec = flatten<pair<double,int>>(T_ec_new);
    for (int i=1;i<T_csr.second.size();i++) { T_rl[i] += T_rl[i-1]; }
    Csr_mtx T_csr_new = make_pair(T_ec,T_rl);
    cout << "leaving regularise_alt()..." << endl;
    return T_csr_new;
}

