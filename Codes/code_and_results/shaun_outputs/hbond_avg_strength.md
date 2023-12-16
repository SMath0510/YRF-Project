# Bond Lifetime Analysis

## Introduction

This document outlines a Python script for analyzing the lifetime of bonds in a molecular dynamics simulation. The script uses NetworkX for graph operations and calculates the average lifetime of each bond over multiple frames.

## Variables and Initialization

```python
old_hbond_list = []
latest_hbond_list = []
average_hbond_life_mapping = {}
```

- `old_hbond_list`: List to store the previous frame's hydrogen bonds.
- `latest_hbond_list`: List to store the current frame's hydrogen bonds.
- `average_hbond_life_mapping`: Dictionary to store the start and end frames of each bond, along with a list of their lifetimes.

## `update_mapping` Function

```python
def update_mapping(frame_no, num_nodes, max_frames=2):
    # Update bond information for the current frame

    for key in latest_hbond_list:
        # Check if the bond is new or if its start frame is not recorded
        if key not in average_hbond_life_mapping.keys() or average_hbond_life_mapping[key]['Start'] is None:
            average_hbond_life_mapping[key] = {'Start': frame_no, 'End': None, 'Time_List': []}

        # Update the end frame for the current frame
        average_hbond_life_mapping[key]['End'] = frame_no

    # Update bond information for the previous frame
    for key in old_hbond_list:
        if key not in average_hbond_life_mapping.keys():
            average_hbond_life_mapping[key] = {'Start': None, 'End': None, 'Time_List': []}

        # Check if the bond is no longer present in the latest frame
        if key not in latest_hbond_list:
            # Update the end frame and calculate the lifetime
            average_hbond_life_mapping[key]['End'] = frame_no
            average_hbond_life_mapping[key]['Time_List'].append(
                average_hbond_life_mapping[key]['End'] - average_hbond_life_mapping[key]['Start']
            )

    # If it's the last frame, update end frames for all bonds
    if (frame_no + 1) == max_frames:
        for key in latest_hbond_list:
            average_hbond_life_mapping[key]['End'] = frame_no + 1
            average_hbond_life_mapping[key]['Time_List'].append(
                average_hbond_life_mapping[key]['End'] - average_hbond_life_mapping[key]['Start']
            )
```

- The function `update_mapping` updates the bond information for the current and previous frames.
- For each bond, it checks whether it's new, still present, or no longer present.
- If it's the last frame, it updates the end frames for all bonds.

## Updating Bond Lists

```python
old_hbond_list = copy.deepcopy(latest_hbond_list)
latest_hbond_list = copy.deepcopy(hbnet.edges)
```

- Save copies of the current frame's bonds as the old bonds and update the latest bonds.

## Calculating Average Bond Lifetime and Saving to CSV

```python
for key in average_hbond_life_mapping.keys():
    len_ = len(average_hbond_life_mapping[key]['Time_List'])
    if len_ == 0:
        continue
    sum_ = sum(average_hbond_life_mapping[key]['Time_List'])
    average = sum_ / len_
    
    file9.write(f'{key}, {average}\n')
```

- Iterate over the recorded bonds in `average_hbond_life_mapping`.
- Calculate the average lifetime for each bond based on the recorded time intervals.
- Save the results to a CSV file.