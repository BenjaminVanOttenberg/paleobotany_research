import numpy as np
import matplotlib.pyplot as plt

# Initialize generation 0
def prometheus_generation(n):
    gen0 = np.array([[i % 2, -1, -1, i, i] for i in range(n)])  # Each individual is its own ancestor
    return gen0

# Create the next generation and track descendants
def create_next_generation(current_gen, n, descendant_counter):
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
        
        # Increment descendant count for original ancestors
        descendant_counter[individual[3]] += 1  # Update count for mother's female ancestor
        descendant_counter[individual[4]] += 1  # Update count for father's male ancestor
        
        new_gen.append(individual)
    
    return np.array(new_gen)

# Simulate multiple generations and track descendant counts
def simulate_generations(initial_gen_size, num_generations, gen_size_per_step):
    gen0 = prometheus_generation(initial_gen_size)
    current_gen = gen0
    num_ancestors = len(gen0)
    
    # Track descendants over generations
    descendant_history = np.zeros((num_generations, num_ancestors))
    
    for gen in range(num_generations):
        descendant_counter = np.zeros(num_ancestors)
        current_gen = create_next_generation(current_gen, gen_size_per_step, descendant_counter)
        descendant_history[gen] = descendant_counter  # Record the descendant counts for this generation
    
    return gen0, descendant_history

# Plot the descendant history over generations
def plot_descendant_history(descendant_history):
    generations = np.arange(descendant_history.shape[0])
    
    plt.figure(figsize=(12, 6))
    plt.stackplot(generations, descendant_history.T, labels=[f'Ancestor {i}' for i in range(descendant_history.shape[1])])
    plt.xlabel('Generations')
    plt.ylabel('Number of Living Descendants')
    plt.title('Living Descendants of Each Ancestor Over Generations')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.tight_layout()
    plt.show()

# Parameters
initial_gen_size = 400       # Initial ancestor generation size
num_generations =500        # Number of generations to simulate
gen_size_per_step = 400      # Number of individuals per new generation

# Run simulation and plot results
gen0, descendant_history = simulate_generations(initial_gen_size, num_generations, gen_size_per_step)
plot_descendant_history(descendant_history)

