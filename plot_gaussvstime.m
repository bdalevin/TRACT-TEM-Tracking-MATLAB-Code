function plot_gaussvstime(TrackedJumpsGauss, ID, Scale)

JumpsID = TrackedJumpsGauss(:,:,ID);
Average_X_Position =  mean(JumpsID(:,2))
Average_Y_Position =  mean(JumpsID(:,4))
t = [0.1:0.1:1];
plot(t,JumpsID(:,10)*1000*Scale,'-o');
xlabel('Time (s)');
ylabel('Jump Distance (pm)');
end