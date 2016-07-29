Vid1Start = 31; % If 31, it will start at frame 1 without clipping
Vid1End = 772; % beginning first movement artefact
Vid2Start = 796; % end movement artefact
% Vid2End = 2037; % beginning first movement artefact
% Vid3Start = 2075;
% Vid3End = 3665;
% Vid4Start = 3706;
% Vid4End = 2812;
% Vid5Start = 2853;


% video1 = video.frames(Vid1Start-30:Vid1End+30);
% % video2 = video.frames(Vid2Start-30:Vid2End+30);
% % video3 = video.frames(Vid3Start-30:Vid3End+30);
% % video4 = video.frames(Vid4Start-30:Vid4End+30);
% video2 = video.frames(Vid2Start-30:end);

video1 = video_reg(:,:,1:722+30);
video2 = video_reg(:,:,796-30:3600);

save('SanneSleep28_lny19rb_2min_iso(13)part1','video1')
save('SanneSleep28_lny19rb_2min_iso(13)part2','video2')
% save('LNY19RB_RH_9-16-15_SanneSleep8(6)part3','video3')
% save('LNY19RB_RH_9-16-15_SanneSleep8(6)part4','video4')
% save('SanneSleep19_LBr40_112815_(4)part5','video5')