import re

location = 'complement(join(4686458..4687873,>1..>60, 76543))'

raw_entries = re.split(',\s?', location) # dont even trust them to comma-separate without spaces
print(raw_entries)
for entry in raw_entries:
    raw_loc = re.search(r'[0-9]+(\.\.[<>]?[0-9]*)?', entry).group()
    loc_coords = re.split(r'\.\.[<>]?', raw_loc)
    if len(loc_coords) == 1:
        loc_coords.append(loc_coords[0])
    elif len(loc_coords) == 0:
        print('something went wrong')
        print(location)
    print(loc_coords)