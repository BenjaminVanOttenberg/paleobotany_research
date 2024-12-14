import numpy as np

# Individual schema:
# [0] gender: 0 (female), 1 (male)
# [1] mother index
# [2] father index
# [3] female ancestor index from gen0
# [4] male ancestor index from gen0

def prometheus_generation(n):
    """
    Initialize generation 0 with `n` individuals where each individual has no parents 
    and no ancestors set (these will be the ancestors).
    """
    gen0 = np.array([[i % 2, -1, -1, i, i] for i in range(n)])  # Each individual is its own ancestor
    return gen0

def create_next_generation(current_gen, n):
    """
    Create the next generation `n` individuals, linking them to parents in the previous generation.
    Each individual will inherit female ancestors from the mother and male ancestors from the father.
    """
    new_gen = []
    
    for _ in range(n):
        # Randomly select a mother and a father
        mother_idx = np.random.choice(np.where(current_gen[:, 0] == 0)[0])
        father_idx = np.random.choice(np.where(current_gen[:, 0] == 1)[0])
        
        # Set gender randomly
        gender = np.random.randint(2)
        
        # Link parents and ancestors
        individual = [
            gender,                     # Gender
            mother_idx,                 # Mother index
            father_idx,                 # Father index
            current_gen[mother_idx, 3], # Female ancestor from mother
            current_gen[father_idx, 4]  # Male ancestor from father
        ]
        
        new_gen.append(individual)
    
    return np.array(new_gen)

# Example usage
gen0 = prometheus_generation(100)    # Generation 0
gen1 = create_next_generation(gen0, 50)  # Generation 1 with 50 individuals

print("Generation 0 (Ancestors):")
print(gen0)
print("\nGeneration 1 (Offspring with linked ancestors):")
print(gen1)