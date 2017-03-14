import pytraj as pt
from pytraj import datafiles

@profile
def process_traj(traj, bnd_list, chunk_sz=60000):
    chunk_id = 0 
    print(len(bnd_list))
    for chunk in traj.iterchunk(chunk_sz, start=0, stop=-1):            
        chunk_end =  (chunk_id + 1) * chunk_sz if (chunk_id + 1) * chunk_sz < traj.n_frames else traj.n_frames       
        print("block stast, end: ", chunk_id * chunk_sz, chunk_end)
        if chunk_end - chunk_id * chunk_sz > 0:
            bonds_val = np.reshape(np.zeros(chunk_sz * len(bnd_list)), (32, 60000)) 
            #bonds_val = pt.distance(chunk, bnd_list, dtype='ndarray')
            print(bonds_val.shape)
            dummy()
            del bonds_val
        chunk_id += 1
        gc.collect()
