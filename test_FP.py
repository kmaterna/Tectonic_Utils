

import second_focal_plane

np1=[327, 6, 105];
# np1=[132,84,88];
np1=[38, 51, 44];

np2 = second_focal_plane.find_aux_plane(np1[0], np1[1], np1[2]);
print(np2);

