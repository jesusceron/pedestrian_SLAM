% Map shown in the main UI
% figure
grid on

xlim([0,20])
ylim([0, 12])
hold on
xlabel('x(meters)')
ylabel('y(meters)')

% Plot the simulated apartment
rectangle('Position',[5.5 2 12.6 6.6])

% %stairs
rectangle('Position',[2.2 2.2500 1.5 0.3125])
rectangle('Position',[2.2 2.5625 1.5 0.3125])
rectangle('Position',[2.2 2.8750 1.5 0.3125])
rectangle('Position',[2.2 3.1875 1.5 0.3125])
rectangle('Position',[2.2 3.5000 1.5 0.3125])
rectangle('Position',[2.2 3.8125 1.5 0.3125])
rectangle('Position',[2.2 4.1250 1.5 0.3125])
rectangle('Position',[2.2 4.4375 1.5 0.3125])
rectangle('Position',[2.2 4.7500 1.5 0.3125])
rectangle('Position',[2.2 5.0625 1.5 0.3125])
rectangle('Position',[2.2 5.3750 1.5 0.3125])
rectangle('Position',[2.2 5.6875 1.5 0.3125])
rectangle('Position',[2.2 6 1.5 0.3125])
rectangle('Position',[2.2 6.3125 1.5 0.3125])
rectangle('Position',[2.2 6.6250 1.5 0.3125])

%zones without use:
rectangle('Position',[5.5 2 4.3 3], 'FaceColor','k')
rectangle('Position',[5.5 6.4 1.5 2.2], 'FaceColor','k')

% Room and living room
rectangle('Position',[11.8 2 0.8 1.8], 'FaceColor','k') % dining table
rectangle('Position',[13.9 2 0.6 1.4], 'FaceColor','k') % wall btw dining and living room
rectangle('Position',[15.1 2 1.4 0.6], 'FaceColor','k') % Living table
rectangle('Position',[13.9 4.6 4.2 0.4], 'FaceColor','k') % wall btw room and living room
rectangle('Position',[16.5 6 1.6 1.6], 'FaceColor','k') % bed

% Kitchen
rectangle('Position',[8.1 6 3.1 0.6], 'FaceColor','k') % Wall btw bathroom-kitchen and main corridor
rectangle('Position',[12.5 6 1.4 0.6], 'FaceColor','k') % Wall btw kitchen and corridor
rectangle('Position',[9.5 6.6 0.8 2], 'FaceColor','k') % Wall btw bathroom and kitchen
rectangle('Position',[13.1 6.6 0.8 2], 'FaceColor','k') % Wall btw kitchen and room
rectangle('Position',[10.3 8 2.8 0.6], 'FaceColor','k') % Tables of the kitchen