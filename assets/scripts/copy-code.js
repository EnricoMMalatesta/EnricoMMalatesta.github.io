document.addEventListener("DOMContentLoaded", function () {
  const codeBlocks = document.querySelectorAll(".franklin-content pre");

  codeBlocks.forEach(function (pre) {
    const code = pre.querySelector("code");
    if (!code) return;

    // Avoid adding the button twice if the script is re-run.
    if (pre.parentElement.classList.contains("code-block-wrapper")) return;

    const wrapper = document.createElement("div");
    wrapper.className = "code-block-wrapper";

    pre.parentNode.insertBefore(wrapper, pre);
    wrapper.appendChild(pre);

    const button = document.createElement("button");
    button.className = "copy-code-button";
    button.type = "button";
    button.innerText = "Copy";

    button.addEventListener("click", async function () {
      const text = code.innerText;

      try {
        await navigator.clipboard.writeText(text);

        button.innerText = "Copied!";
        setTimeout(function () {
          button.innerText = "Copy";
        }, 1500);
      } catch (err) {
        button.innerText = "Failed";
        setTimeout(function () {
          button.innerText = "Copy";
        }, 1500);
      }
    });

    wrapper.appendChild(button);
  });
});
