function borda_results = borda_ranking(all_permutation_ranks, dataset_names)
% BORDA_RANKING - Calculates a consensus score and rank from multiple ranking results.
%
% Syntax:
%   borda_results = borda_ranking(all_permutation_ranks, dataset_names)
%
% Description:
%   This function uses a normalized Borda count method to determine a stable, overall rank from the results of multiple hierarchy permutations.
%
% Workflow:
%   1.  Point Assignment: For each permutation, awards Borda points based on rank (Rank 1 = N-1 points, Rank N = 0 points).
%   2.  Aggregation: Sums the points for each dataset across all permutations.
%   3.  Normalization: Normalizes the total score to a 0-100% scale, where 100% is a unanimous first place.
%   4.  Final Ranking: Sorts datasets by their normalized score to determine the final consensus rank.
%   5.  Distribution Analysis: Calculates the frequency distribution of ranks (e.g., how often a dataset was Rank 1, 2, etc.) for visualization.
%
% Inputs:
%   all_permutation_ranks - Matrix [num_datasets, num_permutations] with the ranks from each analysis.
%   dataset_names         - Cell array with the names of the datasets.
%
% Outputs:
%   borda_results         - Struct with the consensus score (`.score`), final rank (`.rank`), and the detailed rank distribution (`.rank_distribution`).
%
% Author:   Lukas von Erdmannsdorff
% Date:     12.10.2025
% Version:  1

    [num_datasets, num_perms] = size(all_permutation_ranks);
    
    % Point assignment (Borda Count). Rank 1 -> N-1 points, ..., Rank N -> 0 points.
    points_per_rank = (num_datasets - 1):-1:0;
    
    % Calculate total points per dataset
    total_points = zeros(num_datasets, 1);
    for d = 1:num_datasets
        for p = 1:num_perms
            rank_of_d = all_permutation_ranks(d, p);
            total_points(d) = total_points(d) + points_per_rank(rank_of_d);
        end
    end
    
    % Maximum achievable score (for normalization)
    max_possible_points = (num_datasets - 1) * num_perms;
    
    % Calculate the final "consensus score" as a percentage
    borda_score = (total_points / max_possible_points) * 100;
    
    % Final ranking based on the score
    [~, sorted_indices] = sort(borda_score, 'descend');
    final_borda_rank = zeros(num_datasets, 1);
    final_borda_rank(sorted_indices) = 1:num_datasets;

    % Calculate rank distribution for the lollipop plot
    % Stores how often each dataset achieved each rank
    rank_distribution = cell(num_datasets, 1);
    for d = 1:num_datasets
        all_ranks_for_dataset = all_permutation_ranks(d, :);
        unique_ranks = unique(all_ranks_for_dataset);
        
        counts = zeros(size(unique_ranks));
        for k = 1:numel(unique_ranks)
            counts(k) = sum(all_ranks_for_dataset == unique_ranks(k));
        end
        % Ensure both vectors are concatenated as column vectors
        rank_distribution{d} = [unique_ranks(:), counts(:)];
    end
    
    % Bundle results
    borda_results = struct();
    borda_results.score = borda_score;
    borda_results.rank = final_borda_rank;
    borda_results.rank_distribution = rank_distribution;
    borda_results.dataset_names = dataset_names;
end