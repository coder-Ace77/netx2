/**
 * \ingroup structural
 * \function igraph_betweenness_cutoff
 * \brief Range-limited betweenness centrality.
 *
 * This function computes a range-limited version of betweenness centrality
 * by considering only those shortest paths whose length is no greater
 * then the given cutoff value.
 *
 * \param graph The graph object.
 * \param res The result of the computation, a vector containing the
 *        range-limited betweenness scores for the specified vertices.
 * \param vids The vertices for which the range-limited betweenness centrality
 *        scores will be computed.
 * \param directed Logical, if true directed paths will be considered
 *        for directed graphs. It is ignored for undirected graphs.
 * \param weights An optional vector containing edge weights for
 *        calculating weighted betweenness. No edge weight may be NaN.
 *        Supply a null pointer here for unweighted betweenness.
 * \param cutoff The maximal length of paths that will be considered.
 *        If negative, the exact betweenness will be calculated, and
 *        there will be no upper limit on path lengths.
 * \return Error code:
 *        \c IGRAPH_ENOMEM, not enough memory for
 *        temporary data.
 *        \c IGRAPH_EINVVID, invalid vertex ID passed in
 *        \p vids.
 *
 * Time complexity: O(|V||E|),
 * |V| and
 * |E| are the number of vertices and
 * edges in the graph.
 * Note that the time complexity is independent of the number of
 * vertices for which the score is calculated.
 *
 * \sa \ref igraph_betweenness() to calculate the exact betweenness and
 * \ref igraph_edge_betweenness_cutoff() to calculate the range-limited
 * edge betweenness.
 */
igraph_error_t igraph_stress_cutoff(const igraph_t *graph, igraph_vector_t *res,
                              const igraph_vs_t vids, igraph_bool_t directed,
                              const igraph_vector_t *weights, igraph_real_t cutoff) {

    igraph_integer_t no_of_nodes = igraph_vcount(graph);
    igraph_integer_t no_of_edges = igraph_ecount(graph);
    igraph_adjlist_t adjlist, parents;
    igraph_inclist_t inclist;
    igraph_integer_t source, j, neighbor;
    igraph_stack_int_t S;
    igraph_neimode_t mode = directed ? IGRAPH_OUT : IGRAPH_ALL;
    igraph_vector_t dist;
    /* Note: nrgeo holds the number of shortest paths, which may be very large in some cases,
     * e.g. in a grid graph. If using an integer type, this results in overflow.
     * With a 'long long int', overflow already affects the result for a grid as small as 36*36.
     * Therefore, we use a 'igraph_real_t' instead. While a 'igraph_real_t' holds fewer digits than a
     * 'long long int', i.e. its precision is lower, it is effectively immune to overflow.
     * The impact on the precision of the final result is negligible. The max betweenness
     * is correct to 14 decimal digits, i.e. the precision limit of 'igraph_real_t', even
     * for a 101*101 grid graph. */
    igraph_real_t *nrgeo = 0;
    igraph_real_t *tmpscore;
    igraph_vector_t v_tmpres, *tmpres = &v_tmpres;
    igraph_vit_t vit;

    IGRAPH_CHECK(igraph_i_betweenness_check_weights(weights, no_of_edges));

    if (weights) {
        IGRAPH_CHECK(igraph_inclist_init(graph, &inclist, mode, IGRAPH_NO_LOOPS));
        IGRAPH_FINALLY(igraph_inclist_destroy, &inclist);
    } else {
        IGRAPH_CHECK(igraph_adjlist_init(graph, &adjlist, mode, IGRAPH_NO_LOOPS, IGRAPH_MULTIPLE));
        IGRAPH_FINALLY(igraph_adjlist_destroy, &adjlist);
    }

    IGRAPH_CHECK(igraph_adjlist_init_empty(&parents, no_of_nodes));
    IGRAPH_FINALLY(igraph_adjlist_destroy, &parents);

    IGRAPH_CHECK(igraph_stack_int_init(&S, no_of_nodes));
    IGRAPH_FINALLY(igraph_stack_int_destroy, &S);

    IGRAPH_VECTOR_INIT_FINALLY(&dist, no_of_nodes);

    nrgeo = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(nrgeo, "Insufficient memory for betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, nrgeo);

    tmpscore = IGRAPH_CALLOC(no_of_nodes, igraph_real_t);
    IGRAPH_CHECK_OOM(tmpscore, "Insufficient memory for betweenness calculation.");
    IGRAPH_FINALLY(igraph_free, tmpscore);

    if (igraph_vs_is_all(&vids)) {
        /* result covers all vertices */
        IGRAPH_CHECK(igraph_vector_resize(res, no_of_nodes));
        igraph_vector_null(res);
        tmpres = res;
    } else {
        /* result needed only for a subset of the vertices */
        IGRAPH_VECTOR_INIT_FINALLY(tmpres, no_of_nodes);
    }

    for (source = 0; source < no_of_nodes; source++) {

        /* Loop invariant that is valid at this point:
         *
         * - the stack S is empty
         * - the 'dist' vector contains zeros only
         * - the 'nrgeo' array contains zeros only
         * - the 'tmpscore' array contains zeros only
         * - the 'parents' adjacency list contains empty vectors only
         */

        IGRAPH_PROGRESS("Betweenness centrality: ", 100.0 * source / no_of_nodes, 0);
        IGRAPH_ALLOW_INTERRUPTION();

        /* Conduct a single-source shortest path search from the source node */
        if (weights) {
            IGRAPH_CHECK(igraph_i_sspf_weighted(graph, source, &dist, nrgeo, weights, &S, &parents, &inclist, cutoff));
        } else {
            IGRAPH_CHECK(igraph_i_sspf(source, &dist, nrgeo, &S, &parents, &adjlist, cutoff));
        }

        /* Aggregate betweenness scores for the nodes we have reached in this
         * traversal */
        while (!igraph_stack_int_empty(&S)) {
            igraph_integer_t actnode = igraph_stack_int_pop(&S);
            igraph_vector_int_t *neis = igraph_adjlist_get(&parents, actnode);
            igraph_integer_t nneis = igraph_vector_int_size(neis);
            // igraph_real_t coeff = (1 + tmpscore[actnode]) / nrgeo[actnode];

            // for (j = 0; j < nneis; j++) {
            //     neighbor = VECTOR(*neis)[j];
            //     tmpscore[neighbor] += nrgeo[neighbor] * coeff;
            // }

            // if (actnode != source) {
            //     VECTOR(*tmpres)[actnode] += tmpscore[actnode];
            // }
            if (actnode != source) {
                VECTOR(*tmpres)[actnode] += nrgeo[actnode];
            }

            /* Reset variables to ensure that the 'for' loop invariant will
             * still be valid in the next iteration */

            VECTOR(dist)[actnode] = 0;
            nrgeo[actnode] = 0;
            // tmpscore[actnode] = 0;
            igraph_vector_int_clear(neis);
        }

    } /* for source < no_of_nodes */

    /* Keep only the requested vertices */
    if (!igraph_vs_is_all(&vids)) {
        IGRAPH_CHECK(igraph_vit_create(graph, vids, &vit));
        IGRAPH_FINALLY(igraph_vit_destroy, &vit);
        IGRAPH_CHECK(igraph_vector_resize(res, IGRAPH_VIT_SIZE(vit)));

        for (j = 0, IGRAPH_VIT_RESET(vit); !IGRAPH_VIT_END(vit);
             IGRAPH_VIT_NEXT(vit), j++) {
            igraph_integer_t node = IGRAPH_VIT_GET(vit);
            VECTOR(*res)[j] = VECTOR(*tmpres)[node];
        }

        igraph_vit_destroy(&vit);
        igraph_vector_destroy(tmpres);
        IGRAPH_FINALLY_CLEAN(2);
    }

    if (!directed || !igraph_is_directed(graph)) {
        igraph_vector_scale(res, 0.5);
    }

    IGRAPH_PROGRESS("Betweenness centrality: ", 100.0, 0);

    IGRAPH_FREE(nrgeo);
    IGRAPH_FREE(tmpscore);
    igraph_vector_destroy(&dist);
    igraph_stack_int_destroy(&S);
    igraph_adjlist_destroy(&parents);
    if (weights) {
        igraph_inclist_destroy(&inclist);
    } else {
        igraph_adjlist_destroy(&adjlist);
    }
    IGRAPH_FINALLY_CLEAN(6);

    return IGRAPH_SUCCESS;
}
