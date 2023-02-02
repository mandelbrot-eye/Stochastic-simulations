# D.Gillespie paper on the exact stochastic simulations algorithm

from numpy import log, random, zeros, array 

# calculate time interval until the next reaction occurs and choose the site
def find_chain(kon, k1, k2, L, t):
    r1 = random.uniform()
    t +=  1/kon*log(1/r1) # at this step a protein connects DNA next time

    pos =  random.randint(1, L+1) # randomly choosen a position on DNA
    is_occupied = False
    
    return [pos, is_occupied, t]

def choose_state(k1, k2, is_occupied):
    if is_occupied:
        k_obst = k1
    else:
        k_obst = k2
        
    return k_obst    


def chanal(a0, L, t):
    r1 = random.uniform()
    t += 1/a0*log(1/r1) 
    react = a0*r1 # chanal of reaction 

    return [react, t]

def attach(kon, k1, k2, L, pos_obst, is_occupied, t):

    k_obst = choose_state(k1, k2, is_occupied)
    a0 = kon + k_obst

    # mimic a kinetic chanal: 
    [react, t] = chanal(a0, L, t)
    
    # go to the different state of the chain as many times as necessary before attach to the chain       
    while (react > kon):
        is_occupied = not is_occupied # change state
        k_obst = choose_state(k1, k2, is_occupied)
        
        a0 = kon + k_obst
        [react, t] = chanal(a0, L, t)

    pos = random.randint(1, L+1) # randomly choosen a position on DNA  

    while (pos == pos_obst and is_occupied):
        pos =  random.randint(1, L+1) # randomly choosen a position on DNA  
    
    return [pos, is_occupied, t]


def bulk_walks(pos, koff, u, kon, k1, k2, L, pos_obst, is_occupied, t):

    k_obst = choose_state(k1, k2, is_occupied)
        
    a0 = 2*u+koff+k_obst
    
    # mimic a kinetic chanal
    [react, t] = chanal(a0, L, t)
    
    # go to the  left        
    if (react <= u):
        pos -= 1 

    # go to the right        
    if (react > u and react <= 2*u):
        pos += 1 
    
    # go to the solution    
    if (react > 2*u and react <= 2*u+koff): # and react <= a0
        [pos, is_occupied, t] = attach(kon, k1, k2, L, pos_obst, is_occupied, t)

    # go to another subchain
    if (react > 2*u+koff and react <= 2*u+koff+k_obst): # and react <= a0
        is_occupied = not is_occupied # change state
        

    return [pos, is_occupied, t]

def bc_n1_walks(pos, koff, u, kon, k1, k2, L, pos_obst, is_occupied, t):

    k_obst = choose_state(k1, k2, is_occupied)
        
    a0 = u+koff+k_obst

    # mimic a kinetic chanal
    [react, t] = chanal(a0, L, t)

    # go to the right        
    if (react <= u):
        pos += 1 
    
    # go to the solution    
    if (react > u and react <= u+koff): # and react <= a0
        [pos, is_occupied, t] = attach(kon, k1, k2, L, pos_obst, is_occupied, t)

    # go to another subchain
    if (react > u+koff and react <= u+koff+k_obst): # and react <= a0
        is_occupied = not is_occupied # change state
        
    return [pos, is_occupied, t]

def bc_nL_walks(pos, koff, u, kon, k1, k2, L, pos_obst, is_occupied, t):

    k_obst = choose_state(k1, k2, is_occupied)

    a0 = u+koff+k_obst # u - go left + koff - detach from the DNA

    # mimic a kinetic chanal
    [react, t] = chanal(a0, L, t)
    
    # go to the left        
    if (react <= u):
        pos -= 1 
    
    # go to the solution    
    if (react > u and react <= u+koff): # and react <= a0
        [pos, is_occupied, t] = attach(kon, k1, k2, L, pos_obst, is_occupied, t)

    # go to another subchain
    if (react > u+koff and react <= u+koff+k_obst): # and react <= a0
        is_occupied = not is_occupied # change state
        

    return [pos, is_occupied, t]


def obst_m1_walks(pos, koff, u, kon, k1, k2, L, pos_obst, is_occupied, t):

    if is_occupied :
        k_obst = k1 # we are on the upper chain
        a0 = u+koff+k_obst 
    else:
        k_obst = k2
        a0 = 2*u+koff+k_obst 

    # mimic a kinetic chanal
    [react, t] = chanal(a0, L, t)

    if is_occupied:
        # go to left on the upper chain        
        if (react <= u):
            pos -= 1
    
        # go to the solution    
        if (react > u and react <= u+koff): # and react <= a0
            [pos, is_occupied, t] = attach(kon, k1, k2, L, pos_obst, is_occupied, t)

        # go to another subchain
        if (react > u+koff and react <= u+koff+k_obst): # and react <= a0
            is_occupied = not is_occupied # change state

    else:
        # go to the left        
        if (react <= u):
            pos -= 1 

        # go to the right        
        if (react > u and react <= 2*u):
            pos += 1 
    
        # go to the solution    
        if (react > 2*u and react <= 2*u+koff): # and react <= a0
            [pos, is_occupied, t] = attach(kon, k1, k2, L, pos_obst, is_occupied, t)

        # go to another subchain
        if (react > 2*u+koff and react <= 2*u+koff+k_obst): # and react <= a0
            is_occupied = not is_occupied # change state
        
    return [pos, is_occupied, t]


def obst_p1_walks(pos, koff, u, kon, k1, k2, L, pos_obst, is_occupied, t):

    if is_occupied :
        k_obst = k1
        a0 = u+koff+k_obst #
    else:
        k_obst = k2
        a0 = 2*u+koff+k_obst 

    # mimic a kinetic chanal
    [react, t] = chanal(a0, L, t)

    if is_occupied:
        # go to left on the upper chain        
        if (react <= u):
            pos += 1
    
        # go to the solution    
        if (react > u and react <= u+koff): # and react <= a0
            [pos, is_occupied, t] = attach(kon, k1, k2, L, pos_obst, is_occupied, t)

        # go to another subchain
        if (react > u+koff and react <= u+koff+k_obst): # and react <= a0
            is_occupied = not is_occupied # change state

    else:
        # go to the left        
        if (react <= u):
            pos -= 1 

        # go to the right        
        if (react > u and react <= 2*u):
            pos += 1 
    
        # go to the solution    
        if (react > 2*u and react <= 2*u+koff): # and react <= a0
            [pos, is_occupied, t] = attach(kon, k1, k2, L, pos_obst, is_occupied, t)

        # go to another subchain
        if (react > 2*u+koff and react <= 2*u+koff+k_obst): # and react <= a0
            is_occupied = not is_occupied # change state
        
    return [pos, is_occupied, t]


def obst_walks(pos, koff, u, kon, k1, k2, L, pos_obst, is_occupied, t):

    a0 = 2*u+koff

    # mimic a kinetic chanal
    [react, t] = chanal(a0, L, t)


    # go to the left        
    if (react <= u):
        pos -= 1 

    # go to the right        
    if (react > u and react <= 2*u):
        pos += 1 
    
    # go to the solution    
    if (react > 2*u and react <= 2*u+koff): # and react <= a0
        [pos, is_occupied, t] = attach(kon, k1, k2, L, pos_obst, is_occupied, t)

        
    return [pos, is_occupied, t]

        
#parameters

L = 1000        # number of sites on DNA
u    = 10.0**5 # [1/s]
#koff = 10.0**(1) # [1/s]
kon  = 0.1*L # [1/s]

k1 = 10**8
k2 = 10

m = L/2
pos_obst = 3*L/4
num_sim = 10

lam = array([10**(-3), 10**(-2), 10**(-1), 1.0, 10.0, 10**2, 10**3, 10**4, 10**5, 10**6, 10**7])
koff = u/lam**2

t0 = zeros(len(lam)) # mfpt
for i in range(len(lam)): # loop over lambda[i]
    t = 0 # start time for new simulation
    for num in range(num_sim): # loop over the number of simulations 

        # protein starts it's search from the solution. 
        [pos, is_occupied, t] = find_chain(kon, k1, k2, L, t)    


        while (pos != m): # protein beyond the target 
            
            # general case (DNA bulk)      
            if (pos != 1 and pos != L and abs(pos-pos_obst) > 1): 
                [pos, is_occupied, t] = bulk_walks(pos, koff[i], u, kon, k1, k2, L, pos_obst, is_occupied, t)

            # b.c n = 1 - first site
            elif (pos == 1 and pos != pos_obst):         
                [pos, is_occupied, t] = bc_n1_walks(pos, koff[i], u, kon, k1, k2, L, pos_obst, is_occupied, t)
                
            # b.c n = L - last site
            elif (pos == L and pos != pos_obst):
                [pos, is_occupied, t] = bc_nL_walks(pos, koff[i], u, kon, k1, k2, L, pos_obst, is_occupied, t)

            # around obstacle
            elif (pos == pos_obst-1):
                [pos, is_occupied, t] = obst_m1_walks(pos, koff[i], u, kon, k1, k2, L, pos_obst, is_occupied, t)

            elif (pos == pos_obst+1):
                [pos, is_occupied, t] = obst_p1_walks(pos, koff[i], u, kon, k1, k2, L, pos_obst, is_occupied, t)

            elif (pos == pos_obst and (not is_occupied)):
                [pos, is_occupied, t] = obst_walks(pos, koff[i], u, kon, k1, k2, L, pos_obst, is_occupied, t)


    t0[i] = t/num_sim
    print ("T0[", i,"] = ", t0[i])

print ("Done")
