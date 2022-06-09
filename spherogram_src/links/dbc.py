def double_branched_cover(L,n=None):
    """
    Input: Plink L where L has one unknot component U 
        such that the geometric intersection number of L-U with U is +/- 2.
        n is the index of the unknont component; last component by default.
    
    Given a (2,2) tangle T, the annular closure of T, together with its axis U, is such a link L. 
    double_branched_cover assumes the axis U is the last component of L.
    
    Output: Double branched cover of the exterior of T.
    """

    M=L.exterior()
    
    num_components=M.num_cusps()
    
    if num_components  < 2:
        raise ValueError("Number of cusps must be at least 2.")
        
    if n==None:
        n=num_components-1
 
    if n not in all_two_component_clasps(L):
        raise Exception("Given unknot component is not a valid clasp.  See all_two_component_clasps(L).")
    
    fill_list=[(2,0)]*(num_components)
    fill_list[n]=(0,0)
    
    
    
    M.dehn_fill(fill_list)
    covers_list=M.covers(2)
    
    #Look for cover with 2 unfilled cusps and all others with filling (1,0)
    good_covers_list=[]
    
    for C in covers_list:
        
        cusp_pairs=C.cusp_info('filling')
        
        #number of (0,0) cusps
        zero_zero_count=0
        #number of (1,0) cusps
        one_zero_count=0
        
        #list of indices where (0,0) cusps occur
        zero_zero_indices=[]
        
        for i in range(len(cusp_pairs)):
            pair=cusp_pairs[i]
            if pair == (1,0):
                one_zero_count+=1
            elif pair == (0,0):
                zero_zero_count+=1
                zero_zero_indices.append(i)
                
            
        if zero_zero_count == 2 and one_zero_count == len(cusp_pairs)-2:
            good_covers_list.append([C,zero_zero_indices])
            
    if len(good_covers_list)>1:
        raise Exception("MULTIPLE GOOD COVERS")
    elif len(good_covers_list)==0:
        raise Exception("NO GOOD COVERS FOUND")
    else:
        
        #cover with two (0,0) cusps and other cusps (1,0)
        good_cover=good_covers_list[0][0]
        
        #location of (0,0) cusps of the good cover
        good_cover_unfilled_indices=good_covers_list[0][1]
        
        #fill one of the unfilled cusps to get double branched cover of the tangle 
        good_cover.dehn_fill((1,0),good_cover_unfilled_indices[0])
        
        #retriangulate as a one-cusped manifold
        good_cover_f=good_cover.filled_triangulation()
        
        return(good_cover_f)

def dbc_tangle(T):
    """
    Input T is a spherogram (2,2)-tangle.
    Output is the double branched cover of the complement of T.
    
    >>>import snappy
    >>>TestTangle = BraidTangle([2, 2, 1, 1, 2, 2, 1, 1, 2]) + IdentityBraid(1)
    >>>TestTangle.annular_closure(-1).exterior().identify()
    [v3383(0,0)(0,0)(0,0), 8^3_1(0,0)(0,0)(0,0), L8a18(0,0)(0,0)(0,0)]
    
    >>>M=dbc_tangle(TestTangle)
    >>>M.identify()
    [m004(0,0), 4_1(0,0), K2_1(0,0), K4a1(0,0), otet02_00001(0,0)]
    """
    #check that T is (2,2)
    if T.boundary != (2,2):
        raise Exception("Not a (2,2)-tangle.")
    
    L=T.annular_closure(-1)
    return double_branched_cover(L)
    
def all_two_component_clasps(link):
    """
    Given a link, search for unknot components U such that L-U crosses U twice as shown ::
    
           |    |
       ,-- | -- | --,
      |    |    |    |
      |    |    |    |
       `------------'
           |    |
           |    |
    Returns list of indices of all such components
    
    # 1 clasp of 2 strands
    >>>L0 = spherogram.Link('DT: [(-18,8,10,16,14,6,20,4),(-2,12)], [1,1,0,1,0,1,0,0,0,1]')

    # 3 clasps of 2 strands
    >>>L1 = spherogram.Link('DT: [(18,14,34,-36,32,44,26,40,30,38,-8,-46,-12,-4,-42,-2,-22),'
                     '(20,-6),(16,-28),(-24,10)],'
                     '[1,0,0,1,1,0,1,1,0,0,1,1,0,0,0,1,1,1,0,0,1,0,1]')
    """
    
    
    acceptable = [[True, True, False, False], [True, False, False, True],
              [False, False, True, True], [False, True, True, False]]
    
    ans = []
    for i, component in enumerate(link.link_components):
        if len(component) == 4:
            over_under = [cs.is_over_crossing() for cs in component]
            if over_under in acceptable:
                ans.append(i)
    return ans
