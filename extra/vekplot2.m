function [pilx,pily]=vekplot2(x,y,u,v,scale,line,mark,fillcolor)
% VEKPLOT2 - plot vectors as arrows 
% [pilx,pily]=vekplot2(x,y,u,v,scale,line,mark,fillcolor)
% 
% This function is essentially the same as quiver, but in VEKPLOT2 the 
% scale is known so that plotting two vector-fields on top of each other
% for comparison, is possible.
%
% scale: the length in the figure when u^2 + v^2 = 1
% pilx,pily : matrices formed so that plot(pilx,pily) gives a vectorplot
%
% See also: QUIVER, QUIVER3, FEATHER, PLOT

% Copyright Per-Olav Rusaas 1997-2001, por@bukharin.hiof.no
%
% Modifications by J. Kristian Sveen, jks@math.uio.no, 15. Nov. 2001

if ischar(x)
    y=load(x);
end 

if nargin==2 | size(x,2)==4
    scale=y; clear y;
    y=x(:,2); u=x(:,3); v=x(:,4); x=x(:,1); 
end

% remove NaN's
% feature added 15-Nov-2001
ii=~isnan(u);
x=x(ii); y=y(ii); u=u(ii); v=v(ii); 
%

x=x(:)'; y=y(:)'; u=u(:)'; v=v(:)';

maksspiss=0.5*scale;  % maks. lengde av spiss-sidene
alfa=pi/6;            % halve vinkelen for spissen

ih=ishold;

x1=x;
x2=x+u*scale;
y1=y;
y2=y+v*scale;
r=sqrt((x2-x1).^2 + (y2-y1).^2);

retcos=[cos(alfa), -sin(alfa) ; sin(alfa), cos(alfa)];

spisslengde=min(0.5*ones(size(r)), 0.25*r);
spisslengde(2,:)=spisslengde(1,:);

lvek1=retcos*[(x1-x2) ; (y1-y2)];
lvek2=retcos'*[(x1-x2) ; (y1-y2)];
lengde=sqrt(lvek1(1,:).^2+lvek1(2,:).^2);
lengde=max(lengde,ones(size(lengde))*1.0e-200);
lengde(2,:)=lengde(1,:);
lvek1=lvek1./lengde .* spisslengde;
lvek2=lvek2./lengde .* spisslengde;

pilx=[x1; x2; x2+lvek1(1,:); x2; x2+lvek2(1,:)];
pily=[y1; y2; y2+lvek1(2,:); y2; y2+lvek2(2,:)];

vecx=[x1; x2; repmat(NaN,size(x1))];
vecy=[y1; y2; repmat(NaN,size(y1))];
px=[x2+lvek1(1,:); x2; x2+lvek2(1,:); repmat(NaN,size(x1))];
py=[y2+lvek1(2,:); y2; y2+lvek2(2,:); repmat(NaN,size(y1))];

if nargin<6 | isempty(line)
    plot(vecx(:),vecy(:),'b-')
    hold on
    plot(px(:),py(:),'b-')
    %plot(pilx,pily,'blue-')
else
    plot(vecx(:),vecy(:),line)
    hold on
    plot(px(:),py(:),line)
    %plot(pilx,pily,line)
end

if nargin==7 & ~isempty(mark)
    hold on
    plot(x,y,mark);
end

if nargin==8 & ~isempty(fillcolor)
    hold on
    ha=plot(x,y,mark);
    set(ha,'markerFaceColor',fillcolor)
end

if ih
    hold on
else
    hold off
end

if nargout==0
    clear pilx pily
end


