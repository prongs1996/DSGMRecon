import Levenshtein
from multiprocessing import Pool
import time
import argparse

def GreedyMedian(cluster):
    return Levenshtein.median(cluster)

def DoubleSidedGreedyMedian(cluster):
    rev = list(map(lambda x: x[len(x)-1:len(x)//2-overlap:-1], cluster))
    cluster = list(map(lambda x: x[:len(x)//2+1+overlap],cluster))
    LtoR = Levenshtein.median(cluster)
    RtoL = Levenshtein.median(rev)
    if len(LtoR) + len(RtoL) < strand_length:
        diff = strand_length - len(LtoR) - len(RtoL)
        # return LtoR + LtoR[-1]*diff + RtoL
        return LtoR + LtoR[-1]*diff + RtoL[::-1]
    return LtoR + RtoL[(strand_length-len(LtoR))-1::-1]

def DoubleSidedBeamGreedyMedian(cluster):
    rev = list(map(lambda x: x[::-1], cluster))
    LtoR = Levenshtein.beamgreedy(cluster,strand_length,b=beam_size)
    RtoL = Levenshtein.beamgreedy(rev,strand_length,b=beam_size)
    return LtoR[:(strand_length+1)//2] + RtoL[(strand_length//2)-1::-1]

def DoubleSidedGreedyMedianRefine(cluster):
    min_length = 1
    rev = list(map(lambda x: x[len(x)-1:len(x)//2-overlap:-1], cluster))
    forward = list(map(lambda x: x[:len(x)//2+1+overlap],cluster))
    LtoR = Levenshtein.median(forward)
    RtoL = Levenshtein.median(rev)
    if len(LtoR) + len(RtoL) < strand_length:
        diff = strand_length - len(LtoR) - len(RtoL)
        # reconstructed = LtoR + LtoR[-1]*diff + RtoL
        reconstructed = LtoR + LtoR[-1]*diff + RtoL[::-1]
    else:
        reconstructed = LtoR + RtoL[(strand_length-len(LtoR))-1::-1]
    curr_run_len = 1
    curr_run_start = 0
    all_long_runs = {}
    for i in range(1,len(reconstructed)):
        if reconstructed[i-1] == reconstructed[i]:
            curr_run_len += 1
        elif curr_run_len >= min_length:
            all_long_runs[curr_run_start] = i
            curr_run_len = 1
            curr_run_start = i
        else:
            curr_run_len = 1
            curr_run_start = i
    if curr_run_len >= min_length:
        all_long_runs[curr_run_start] = i+1
    best_candidate = reconstructed
    curr_total = compute_total_edit_dist(cluster, reconstructed)
    for i in all_long_runs:
        #try insertion
        inserted_strand = reconstructed[:all_long_runs[i]]+reconstructed[i]+reconstructed[all_long_runs[i]:]
        for j in range(len(inserted_strand)):
            if i <= j <= all_long_runs[i]:
                continue
            elif j>0 and (inserted_strand[j] == inserted_strand[j-1]):
                continue
            candidate = inserted_strand[:j]+inserted_strand[j+1:]
            candidate_total = compute_total_edit_dist(cluster, candidate)
            if candidate_total < curr_total:
                curr_total = candidate_total
                best_candidate = candidate
        #try deletion
        deleted_strand = reconstructed[:all_long_runs[i]-1] + reconstructed[all_long_runs[i]:]
        for j in range(len(deleted_strand)):
            if i <= j <= all_long_runs[i]-2:
                continue
            elif j>0 and (deleted_strand[j] == deleted_strand[j-1]):
                continue
            candidate = deleted_strand[:j]+deleted_strand[j]+deleted_strand[j:]
            candidate_total = compute_total_edit_dist(cluster, candidate)
            if candidate_total < curr_total:
                curr_total = candidate_total
                best_candidate = candidate
    return best_candidate

def compute_total_edit_dist(cluster, reconstructed):
    total = 0
    for s in cluster:
        total += Levenshtein.distance(s, reconstructed)
    return total

overlap = 2



### Cluster Parameters
min_coverage = 1
max_coverage = None
strand_length = 110



### Main
if __name__ == '__main__':
    parser = argparse.ArgumentParser('DNA consensus')
    parser.add_argument('--i', type=str, default="our_nanopore_UnderlyingClusters.txt", help="input file")
    parser.add_argument('--r', type=str, default="",help='reference file')
    parser.add_argument('--o', type=str, default="ReconstructedStrands.txt", help="output file")
    parser.add_argument('--ALG', type=int, default=0, help="algorithm number")
    parser.add_argument('--b', type=int, default=10, help='beam size')
    args = parser.parse_args()
    f = [DoubleSidedGreedyMedian, DoubleSidedBeamGreedyMedian, DoubleSidedGreedyMedianRefine][args.ALG]
    beam_size = args.b
    input_file = open("data/"+args.i,'r')
    print(f.__name__+':', args.i)
    data = input_file.readlines()
    clusters = []
    cluster_index = 0
    curr_cluster = []
    skip = set()
    for l in data[1:]:
        if l[1] == 'L':
            if len(curr_cluster) >= min_coverage:
                if max_coverage:
                    clusters.append(curr_cluster[:max_coverage])
                else:
                    clusters.append(curr_cluster)
            else:
                skip.add(cluster_index)
            curr_cluster = []
            cluster_index+=1
        else:
            curr_cluster.append(l.strip())
    if len(curr_cluster) >= min_coverage:
        if max_coverage:
            clusters.append(curr_cluster[:max_coverage])
        else:
            clusters.append(curr_cluster)
    else:
        skip.add(cluster_index)

    pool = Pool()
    print("Num clusters:", len(clusters))
    start = time.time()
    reconstructed_strands = pool.map(f, clusters)
    print("Time:\n\t", time.time() - start)

    ### Write reconstructed strands
    # print("Writing Strands")
    output_file = open("data/"+args.o,'w')
    for line in reconstructed_strands:
        output_file.write(line)
        output_file.write('\n')


    ### Calculating Accuracy
    # print("Calculating Accuracy")
    if args.r:
        reference_file = open("data/"+args.r,'r')
        answer = reference_file.readlines()
        answer = list(map(lambda x:x.strip(), answer))
        offset = 0
        total_num_bases = 0
        total_edit_dist = 0
        total_num_strands = 0
        total_errors = 0
        success = 0
        for i in range(len(answer)):
            if i in skip:
                offset+=1
                continue
            else:
                total_num_strands += 1
                total_edit_dist += Levenshtein.distance(reconstructed_strands[i-offset], answer[i])
                total_num_bases += len(answer[i])
                if reconstructed_strands[i-offset] == answer[i]:
                        success+=1
                for j in range(len(answer[i])):
                    
                    if j>= len(reconstructed_strands[i-offset]):
                        total_errors+= len(answer[i]) - len(reconstructed_strands[i-offset])
                        break
                    elif reconstructed_strands[i-offset][j] != answer[i][j]:
                        total_errors+= 1

        print("(Avg Accuracy, Avg Edit Dist):")
        print('\t',1- (total_errors/total_num_bases), total_edit_dist/total_num_strands)
        print("Success rate:")
        print('\t',success/total_num_strands)
        print('\n--------------------------------------------------\n')