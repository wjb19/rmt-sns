function []=gen_image(file)

figure;
im=load(file);

i=find(im(:,1)==2);
cnt=im(i,:);
plot(cnt(:,8),cnt(:,9),'.',"linewidth",2);
hold on;

i=find(im(:,1)==0);
lons=im(i,:);
nl=unique(lons(:,2));
for i=0:length(nl)
j=find(lons(:,2)==i);
pl=lons(j,:);
plot(pl(:,8),pl(:,9),"linewidth",2);
end

i=find(im(:,1)==1);
lats=im(i,:);
nl=unique(lats(:,2));
for i=0:length(nl)
j=find(lats(:,2)==i);
pl=lats(j,:);
plot(pl(:,8),pl(:,9),"linewidth",2);
end





