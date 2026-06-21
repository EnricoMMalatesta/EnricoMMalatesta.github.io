# Comment feature enabled!

@def title = "Comment feature enabled!"
@def authors = "E. Malatesta"
@def published = "4 January 2026"
@def pt_lang = false
@def rss_pubdate = Date(2026, 01, 04)
@def rss = "Comments!"
@def rss_description = """Comments"""

{{ published }} | **{{ authors }}**

@@post-top-nav
[← All blog posts](/pages/blog/)
@@

## Adding comments to your Franklin blog

As Franklin generates a *static* website, comments are not directly handled by it. As suggested by ChatGPT, luckily enough there exists third-party widgets that can be used to add comment sections. I decided to use [utterances](https://utteranc.es/) but similar alternatives are possible (see for example [giscus](https://giscus.app/)). See also [Ricardo Rosa's website](https://github.com/rmsrosa/rmsrosa.github.io?tab=readme-ov-file).

I would also like to thank ChatGPT for helping me set up this feature. 

{{ blogcomments }}

