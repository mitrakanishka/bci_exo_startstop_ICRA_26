function boxplot_style(data, varargin)
if nargin == 1
    b = boxplot(data,'symbol','');
elseif nargin == 2
    b = boxplot(data,'symbol','', 'Positions', varargin{1});
end

box_color = [0.2 0.4 1.0];
median_color = [0 0 0];
mean_color = [0.8 0 0];
dots_color = [0.7,0.7,0.7];

h=findobj(gca,'tag','Box'); % Box
set(h,'Color',box_color,'LineWidth',2,'MarkerFaceColor','b'); 

h=findobj(gca,'tag','Median'); % Median
set(h,'Color',median_color,'LineWidth',1,'ZData',[10 10]);

whisker_color = box_color; 
h=findobj(gca,'tag','Upper Adjacent Value'); % Whiskers 
set(h,'Color',whisker_color,'LineWidth',2);  
h=findobj(gca,'tag','Lower Adjacent Value'); 
set(h,'Color',whisker_color,'LineWidth',2); 
h=findobj(gca,'tag','Upper Whisker'); 
set(h,'Color',whisker_color,'LineWidth',2); 
h=findobj(gca,'tag','Lower Whisker'); 
set(h,'Color',whisker_color,'LineWidth',2); 

mn = mean(data,1);
x_vals = zeros(2,length(mn));

% Box fill
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    x_d = get(h(j),'XData');
    patch(x_d,get(h(j),'YData'),-1*ones(length(x_d),1),box_color,'FaceAlpha',.2);
    x_vals(1,j) = x_d(1);
    x_vals(2,j) = x_d(3);
end
[~,sort_i] = sort(x_vals(1,:),'ascend');
x_vals = x_vals(:,sort_i);

hold on
for j=1:size(x_vals,2)
    plot(x_vals(:,j),[mn(j) mn(j)],'Color',mean_color,'LineWidth',4,'ZData',[20 20]);
end

% Single subject data
n_sbj = size(data,1);
for i = 1:size(data,2)
    % With newer Matlab try to use 'MarkerAlphaColor',0.2
    if nargin == 1
        scatter(i+0.2*(rand(n_sbj,1)-0.5),data(:,i),70,'filled','MarkerFaceColor',dots_color);
    elseif nargin == 2
        scatter(varargin{1}(i)+0.2*(rand(n_sbj,1)-0.5),data(:,i),70,'filled','MarkerFaceColor',dots_color);
    end
end
