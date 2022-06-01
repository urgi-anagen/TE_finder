#include "FragJoin.h"

//------------------------------------------------------------------------
void FragJoin::align_all(std::list<RangePair> &l,std::list<RangePairSet> &finish_list) {

    std::list<RangePair> list_out_dd;
    std::list<RangePair> list_out_cd;
    std::list<RangePair> list_out_dc;
    std::list<RangePair> list_out_cc;

    splitFromStrand(l,list_out_dd,list_out_cd,list_out_dc,list_out_cc);

    align_dd(list_out_dd,finish_list);
    align_cd(list_out_cd,finish_list);
    align_dc(list_out_dc,finish_list);
    align_cc(list_out_cc,finish_list);
}
//------------------------------------------------------------------------
void FragJoin::splitFromStrand(std::list<RangePair>& list_in,
                                    std::list<RangePair>& list_out_dd,
                                    std::list<RangePair>& list_out_cd,
                                    std::list<RangePair>& list_out_dc,
                                    std::list<RangePair>& list_out_cc)
                                    {
    auto it_rp_in=list_in.begin();
    while(it_rp_in!=list_in.end())
    {
        if (it_rp_in->getRangeQ().isPlusStrand() && it_rp_in->getRangeS().isPlusStrand()){
            list_out_dd.push_back(*it_rp_in);
            it_rp_in=list_in.erase(it_rp_in);
            continue;
        }

        if (it_rp_in->getRangeQ().isPlusStrand() && !it_rp_in->getRangeS().isPlusStrand()){
            list_out_dc.push_back(*it_rp_in);
            it_rp_in=list_in.erase(it_rp_in);
            continue;
        }
        if (!it_rp_in->getRangeQ().isPlusStrand() && it_rp_in->getRangeS().isPlusStrand()){
            list_out_cd.push_back(*it_rp_in);
            it_rp_in=list_in.erase(it_rp_in);
            continue;
        }
        if (!it_rp_in->getRangeQ().isPlusStrand() && !it_rp_in->getRangeS().isPlusStrand()){
            list_out_cc.push_back(*it_rp_in);
            it_rp_in=list_in.erase(it_rp_in);
            continue;
        }
        it_rp_in++;
    }
}
//------------------------------------------------------------------------
void FragJoin::join_path2rp_list(std::list<std::list<RangePair>::iterator>& join_path, std::list<RangePair> &l_in, std::list<RangePairSet> &l_out) {
    if (!join_path.empty()) {
        RangePairSet newRange = *(join_path.front());
        std::list<RangePair> a_path;
        for (std::list< std::list<RangePair>::iterator >::iterator r = join_path.begin();
             r != join_path.end(); r++) {
            a_path.push_back(*(*r));
            l_in.erase(*r);
        }
        newRange.setRpsFromRpList(a_path);
        newRange.computeScoreWithDynaProg(mism_pen, gapo_pen, gape_pen);
        l_out.push_back(newRange);
    }
}
//------------------------------------------------------------------------
void FragJoin::align_dd(std::list<RangePair> &l_in, std::list<RangePairSet> &l_out) {
    while(!l_in.empty()){
        std::list<bound> I; // list of horizontal rectangle coordinates (left and right corners)
        std::vector<long long> V; // vector of the best chain scores that finish in a given rectangle
        V.push_back(0);
        unsigned k = 0;
        unsigned nb_frag = 0;
        for (std::list<RangePair>::iterator i = l_in.begin(); i != l_in.end(); i++) {
            nb_frag++;
            bound b;
            unsigned o = over;
            if (i->getRangeQ().getLength() / 2 < over
                || i->getRangeS().getLength() / 2 < over)
                o = std::min(i->getRangeQ().getLength() / 2, i->getRangeS().getLength() / 2);
            b.end_point = i->getRangeQ().getStart() + o;
            b.y_high = i->getRangeS().getStart() + o;
            b.y_low = i->getRangeS().getEnd() - o;
            b.left = true;
            b.iter_range_pair = i;
            b.rect = ++k;
            I.push_back(b);

            b.end_point = i->getRangeQ().getEnd() - o;
            b.y_high = i->getRangeS().getStart() + o;
            b.y_low = i->getRangeS().getEnd() - o;
            b.left = false;
            b.iter_range_pair = i;
            b.rect = k;
            I.push_back(b);

            V.push_back(0);
        }
        if (nb_frag > 0)
        {
            std::list<std::list<RangePair>::iterator> join_path;
            align(nb_frag,I,V,join_path);
            join_path2rp_list(join_path, l_in, l_out);
        }
    }
}

//------------------------------------------------------------------------
void FragJoin::align_dc(std::list<RangePair> &l_in, std::list<RangePairSet> &l_out) {
    while(!l_in.empty()){
        std::list<bound> I; // list of horizontal rectangle coordinates (left and right corners)
        std::vector<long long> V; // vector of the best chain scores that finish in a given rectangle
        V.push_back(0);
        unsigned k = 0;
        unsigned nb_frag = 0;
        unsigned long max_coord = 0;
        for (std::list<RangePair>::iterator i = l_in.begin(); i != l_in.end(); i++) {
            unsigned long m = std::max(i->getRangeS().getStart(), i->getRangeS().getEnd());
            max_coord = m > max_coord ? m : max_coord;
        }
        max_coord++;
        for (std::list<RangePair>::iterator i = l_in.begin(); i != l_in.end(); i++) {
            nb_frag++;
            bound b;
            unsigned o = over;
            if (i->getRangeQ().getLength() / 2 < over
                || i->getRangeS().getLength() / 2 < over)
                o = std::min(i->getRangeQ().getLength() / 2, i->getRangeS().getLength() / 2);
            b.end_point = i->getRangeQ().getStart() + o;
            b.y_high = max_coord - i->getRangeS().getStart() + o;
            b.y_low = max_coord - i->getRangeS().getEnd() - o;
            b.left = true;
            b.iter_range_pair = i;
            b.rect = ++k;
            I.push_back(b);

            b.end_point = i->getRangeQ().getEnd() - o;
            b.y_high = max_coord - i->getRangeS().getStart() + o;
            b.y_low = max_coord - i->getRangeS().getEnd() - o;
            b.left = false;
            b.iter_range_pair = i;
            b.rect = k;
            I.push_back(b);

            V.push_back(0);
        }
        if (nb_frag > 0)
        {
            std::list<std::list<RangePair>::iterator> join_path;
            align(nb_frag, I,V,join_path);
            join_path2rp_list(join_path, l_in, l_out);
        }
    }
}

//------------------------------------------------------------------------
void FragJoin::align_cd(std::list<RangePair> &l_in, std::list<RangePairSet> &l_out) {
    while(!l_in.empty()){
        std::list<bound> I; // list of horizontal rectangle coordinates (left and right corners)
        std::vector<long long> V; // vector of the best chain scores that finish in a given rectangle
        V.push_back(0);
        unsigned k = 0;
        unsigned nb_frag = 0;
        unsigned long max_coord = 0;
        for (std::list<RangePair>::iterator i = l_in.begin(); i != l_in.end(); i++) {
            unsigned long m = std::max(i->getRangeQ().getStart(), i->getRangeQ().getEnd());
            max_coord = m > max_coord ? m : max_coord;
        }
        max_coord++;
        for (std::list<RangePair>::iterator i = l_in.begin(); i != l_in.end(); i++) {
            nb_frag++;
            bound b;
            unsigned o = over;
            if (i->getRangeQ().getLength() / 2 < over
                || i->getRangeS().getLength() / 2 < over)
                o = std::min(i->getRangeQ().getLength() / 2, i->getRangeS().getLength() / 2);
            b.end_point = max_coord - i->getRangeQ().getStart() + o;
            b.y_high = i->getRangeS().getStart() + o;
            b.y_low = i->getRangeS().getEnd() - o;
            b.left = true;
            b.iter_range_pair = i;
            b.rect = ++k;
            I.push_back(b);

            b.end_point = max_coord - i->getRangeQ().getEnd() - o;
            b.y_high = i->getRangeS().getStart() + o;
            b.y_low = i->getRangeS().getEnd() - o;
            b.left = false;
            b.iter_range_pair = i;
            b.rect = k;
            I.push_back(b);

            V.push_back(0);
        }
        if (nb_frag > 0)
        {
            std::list<std::list<RangePair>::iterator> join_path;
            align(nb_frag, I,V,join_path);
            join_path2rp_list(join_path, l_in, l_out);
        }
    }
}

//------------------------------------------------------------------------
void FragJoin::align_cc(std::list<RangePair> &l_in, std::list<RangePairSet> &l_out) {
    while (!l_in.empty()) {
        std::list<bound> I; // list of horizontal rectangle coordinates (left and right corners)
        std::vector<long long> V; // vector of the best chain scores that finish in a given rectangle
        V.push_back(0);
        unsigned k = 0;
        unsigned nb_frag = 0;
        unsigned long max_coordQ = 0;
        for (std::list<RangePair>::iterator i = l_in.begin(); i != l_in.end(); i++) {
            unsigned long m = std::max(i->getRangeQ().getStart(), i->getRangeQ().getEnd());
            max_coordQ = m > max_coordQ ? m : max_coordQ;
        }
        max_coordQ++;
        unsigned long max_coordS = 0;
        for (std::list<RangePair>::iterator i = l_in.begin(); i != l_in.end(); i++) {
            unsigned long m = std::max(i->getRangeQ().getStart(), i->getRangeQ().getEnd());
            max_coordS = m > max_coordS ? m : max_coordS;
        }
        max_coordS++;
        for (std::list<RangePair>::iterator i = l_in.begin(); i != l_in.end(); i++) {
            nb_frag++;
            bound b;
            unsigned o = over;
            if (i->getRangeQ().getLength() / 2 < over
                || i->getRangeS().getLength() / 2 < over)
                o = std::min(i->getRangeQ().getLength() / 2, i->getRangeS().getLength() / 2);
            b.end_point = max_coordQ - i->getRangeQ().getStart() + o;
            b.y_high = max_coordS - i->getRangeS().getStart() + o;
            b.y_low = max_coordS - i->getRangeS().getEnd() - o;
            b.left = true;
            b.iter_range_pair = i;
            b.rect = ++k;
            I.push_back(b);

            b.end_point = max_coordQ - i->getRangeQ().getEnd() - o;
            b.y_high = max_coordS - i->getRangeS().getStart() + o;
            b.y_low = max_coordS - i->getRangeS().getEnd() - o;
            b.left = false;
            b.iter_range_pair = i;
            b.rect = k;
            I.push_back(b);

            V.push_back(0);
        }
        if (nb_frag > 0)
        {
            std::list<std::list<RangePair>::iterator> join_path;
            align(nb_frag, I,V,join_path);
            join_path2rp_list(join_path, l_in, l_out);
        }
    }
}
//------------------------------------------------------------------------
void FragJoin::align(unsigned nb_frag, std::list<bound>& I, std::vector<long long>& V, std::list<std::list<RangePair>::iterator>& join_path) {
    //Gusfield, page 328
    std::vector<std::list<std::list<RangePair>::iterator> > path(nb_frag + 1);
    unsigned end_best_path = 0;
    V[0] = 0;
    std::list<triple> L; //list of vertical rectangle coordinates (kept sorted)
    I.sort(bound::less);
    std::list<bound>::iterator i = I.begin();
    while (i != I.end()) // for each hz rect coordinate
    {
        unsigned k = i->rect;
        if (i->left) //left coordinate: compute the score and chain
        {
            unsigned long hk = i->y_high;

            long long max_score = std::numeric_limits<long long>::min();
            std::list<triple>::iterator j = L.begin();
            std::list<triple>::iterator max_j = L.end();
            while (j != L.end() && j->lower < hk) // modif Hadi: run_test_search_wSW for ALL rectangle above the current
            {
                long long diag1 =
                        (long long) j->iter_range_pair->getRangeQ().getEnd()
                        - (long long) j->iter_range_pair->getRangeS().getEnd();
                long long diag2 =
                        (long long) i->iter_range_pair->getRangeQ().getStart()
                        - (long long) i->iter_range_pair->getRangeS().getStart();
                long long gap = llabs(diag1 - diag2);

                long long mismatch;
                if (diag1 > diag2)
                    mismatch =
                            (long long) i->iter_range_pair->getRangeQ().getStart()
                            - (long long) j->iter_range_pair->getRangeQ().getEnd();
                else
                    mismatch =
                            (long long) i->iter_range_pair->getRangeS().getStart()
                            - (long long) j->iter_range_pair->getRangeS().getEnd();

                mismatch = mismatch > 0 ? mismatch : 0;

                long long score = i->iter_range_pair->getScore()
                                  + V[j->rect]
                                  - (long long) rint(gapo_pen + gape_pen * gap)
                                  - (long long) rint(mismatch * mism_pen);

                if (score > max_score) {
                    max_score = score;
                    max_j = j;
                }
                j++;
            }

            if (max_j == L.end()) {
                V[k] = i->iter_range_pair->getScore();
                path[k].push_back(i->iter_range_pair);
            } else {
                if (max_score >= i->iter_range_pair->getScore()) {
                    V[k] = max_score;
                    path[k].insert(path[k].end(), path[max_j->rect].begin(),
                                   path[max_j->rect].end());
                    path[k].push_back(i->iter_range_pair);
                } else {
                    V[k] = i->iter_range_pair->getScore();
                    path[k].push_back(i->iter_range_pair);
                }
            }
            if (V[k] > V[end_best_path]) end_best_path = k;
        } else // right coordinate: add vertical coord as a triple
        {
            unsigned long lk = i->y_low;
            std::list<triple>::iterator j = std::lower_bound(L.begin(), L.end(), lk,
                                                             triple::lessI);
            if (j != L.end()) j--;
            long long Vj = 0;
            unsigned long lj = 0;
            if (j != L.end()) {
                Vj = V[j->rect];
                lj = j->lower;
            }
            if (lj < lk || (lj == lk && V[k] > Vj)) {
                triple t;
                t.lower = i->y_low;
                t.path_val = V[k];
                t.iter_range_pair = i->iter_range_pair;
                t.rect = k;
                L.insert(std::lower_bound(L.begin(), L.end(), t, triple::less), t);
            }
            std::list<triple>::iterator n = L.begin();
            while (n != L.end())
                if (V[k] > V[n->rect] && n->lower <= lk
                    && (i->iter_range_pair->getRangeQ().getStart()
                        - n->iter_range_pair->getRangeQ().getEnd()
                        > connect_dist_limit
                        ||
                        i->iter_range_pair->getRangeS().getStart()
                        - n->iter_range_pair->getRangeS().getEnd()
                        > connect_dist_limit
                    )
                    // modif Hadi to take into account above modifications
                    //and set limits to run_test_search_wSW for above rectangles
                        )
                    n = L.erase(n);
                else n++;

        }
        i++;
    }
    join_path = path[end_best_path];
};









